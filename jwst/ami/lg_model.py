import logging
import numpy as np

from . import leastsqnrm as leastsqnrm
from . import analyticnrm2
from . import utils
from . import mask_definition_ami

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())
log.setLevel(logging.DEBUG)


m = 1.0
mm = 1.0e-3 * m
um = 1.0e-6 * m
mas = 1.0e-3 / (60 * 60 * 180 / np.pi)  # in radians


class LgModel:
    """
    A class for conveniently dealing with an "NRM object.

    This should be able to take an NRMDefinition object for mask geometry.
    Defines mask geometry and detector-scale parameters.
    Simulates PSF (broadband or monochromatic)
    Builds a fringe model - either by user definition, or automated to data
    Fits model to data by least squares
    Masks: jwst_ami (formerly jwst_g7s6c)
    Algorithm documented in Greenbaum, A. Z., Pueyo, L. P., Sivaramakrishnan,
    A., and Lacour, S., Astrophysical Journal vol. 798, Jan 2015.
    First written by Alexandra Greenbaum in 2014.
    """

    def __init__(
        self,
        nrm_model,
        mask=None,
        holeshape="hex",
        pixscale=None,
        over=1,
        pixweight=None,
        phi=None,
        chooseholes=None,
        affine2d=None,
        **kwargs,
    ):
        """
        Set attributes of LgModel class.

        Parameters
        ----------
        nrm_model : NRMModel datamodel
            Datamodel containing mask geometry information

        mask : str
            Keyword for built-in values

        holeshape : str
           Shape of apertures, default="hex"

        pixscale : float
           Initial estimate of pixel scale in radians

        over : int
           Oversampling factor

        pixweight : 2D float array, default None
            Weighting array

        phi : float 1D array
            Distance of fringe from hole center in units of waves

        chooseholes : list of strings, default None
            E.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask
            If None, use the real seven-hole mask.

        affine2d : Affine2d object
            Affine2d object

        **kwargs : dict
            Keyword arguments
            debug: boolean
                if set, print debug
        """
        if "debug" in kwargs:
            self.debug = kwargs["debug"]
        else:
            self.debug = False

        self.holeshape = holeshape
        self.pixel = pixscale  # det pix in rad (square)
        self.over = over
        self.pixweight = pixweight

        # get these from mask_definition_ami instead
        if mask is None:
            log.info("Using JWST AMI mask geometry from LgModel")
            mask = mask_definition_ami.NRMDefinition(
                nrm_model, maskname="jwst_ami", chooseholes=chooseholes
            )
        elif isinstance(mask, str):
            mask = mask_definition_ami.NRMDefinition(
                nrm_model, maskname=mask, chooseholes=chooseholes
            )  # retain ability to possibly  use other named masks, for now
        self.ctrs = mask.ctrs
        self.d = mask.hdia
        self.D = mask.active_D

        self.N = len(self.ctrs)
        self.fmt = "%10.4e"

        # get closest in time OPD from WebbPSF?

        if phi:  # meters of OPD at central wavelength
            if phi == "perfect":
                self.phi = np.zeros(self.N)  # backwards compatibility
            else:
                self.phi = phi
        else:
            self.phi = np.zeros(self.N)

        self.chooseholes = chooseholes

        # affine2d property not to be changed in LgModel - create a new
        #     instance instead
        # Save affine deformation of pupil object or create a no-deformation
        #     object.
        # We apply this when sampling the PSF, not to the pupil geometry.

        if affine2d is None:
            self.affine2d = utils.Affine2d(
                mx=1.0, my=1.0, sx=0.0, sy=0.0, xo=0.0, yo=0.0, name="Ideal"
            )
        else:
            self.affine2d = affine2d

    def simulate(self, fov=None, bandpass=None, over=None, psf_offset=(0, 0)):
        """
        Simulate a detector-scale psf.

        Use parameters input from the call and
        already stored in the object, and generate a simulation fits header
        storing all of the  parameters used to generate that psf.  If the input
        bandpass is one number it will calculate a monochromatic psf.

        Parameters
        ----------
        fov : int, default=None
            Number of detector pixels on a side

        bandpass : 2D float array, default=None
            Array of the form: [(weight1, wavl1), (weight2, wavl2), ...]

        over : int
            Oversampling factor

        psf_offset : detector pixels
            Center offset from center of array

        Returns
        -------
        Object's 'psf' : float 2D array
            Simulated psf
        """
        # First set up conditions for choosing various parameters
        self.bandpass = bandpass

        if over is None:
            over = 1  # ?  Always comes in as integer.

        self.psf_over = np.zeros((over * fov, over * fov))
        nspec = 0
        # accumulate polychromatic oversampled psf in the object

        for w, l in bandpass:  # w: wavelength's weight, l: lambda (wavelength)
            self.psf_over += w * analyticnrm2.psf(
                self.pixel,  # det pixel, rad
                fov,  # in detpix number
                over,
                self.ctrs,
                self.d,
                l,
                self.phi,
                psf_offset,  # det pixels
                self.affine2d,
                shape=self.holeshape,
            )
            # offset signs fixed to agree w/DS9, +x shifts ctr R, +y shifts up
            nspec += 1

        # store the detector pixel scale psf in the object
        self.psf = utils.rebin(self.psf_over, (over, over))

        return self.psf

    def make_model(self, fov=None, bandpass=None, over=1, psf_offset=(0, 0), pixscale=None):
        """
        Generate the fringe model.

        Use the attributes of the object with a bandpass that is either a single
        wavelength or a list of tuples of the form
        [(weight1, wavl1), (weight2, wavl2),...].  The model is
        a collection of fringe intensities, where nholes = 7 means the model
        has a @D slice for each of 21 cosines, 21 sines, a DC-like, and a flux
        slice for a toal of 44 2D slices.

        Parameters
        ----------
        fov : int, default=None
            Number of detector pixels on a side

        bandpass : 2D float array, default=None
            Array of the form: [(weight1, wavl1), (weight2, wavl2), ...]

        over : int
           Cversampling factor

        psf_offset : detector pixels
            Center offset from center of array

        pixscale : float, default=None
            Pixel scale

        Returns
        -------
        Object's 'model': fringe model
            Generated fringe model
        """
        if fov:
            self.fov = fov

        self.over = over

        if hasattr(self, "pixscale_measured"):
            if self.pixscale_measured is not None:
                self.modelpix = self.pixscale_measured

        if pixscale is None:
            self.modelpix = self.pixel
        else:
            self.modelpix = pixscale

        self.modelctrs = self.ctrs

        # The model shape is (fov) x (fov) x (# solution coefficients)
        # the coefficient refers to the terms in the analytic equation
        # There are N(N-1) independent pistons, double-counted by cosine
        # and sine, one constant term and a DC offset.

        self.model = np.zeros((self.fov, self.fov, self.N * (self.N - 1) + 2))
        self.model_beam = np.zeros((self.over * self.fov, self.over * self.fov))
        self.fringes = np.zeros(
            (self.N * (self.N - 1) + 1, self.over * self.fov, self.over * self.fov)
        )

        for w, l in bandpass:  # w: weight, l: lambda (wavelength)
            # model_array returns the envelope and fringe model as a list of
            #   oversampled fov x fov slices
            pb, ff = analyticnrm2.model_array(
                self.modelctrs,
                l,
                self.over,
                self.modelpix,
                self.fov,
                self.d,
                shape=self.holeshape,
                psf_offset=psf_offset,
                affine2d=self.affine2d,
            )

            log.debug(f"Passed to model_array: psf_offset: {psf_offset}")
            log.debug(f"Primary beam in the model created: {pb}")
            self.model_beam += pb
            self.fringes += ff

            # this routine multiplies the envelope by each fringe "image"
            self.model_over = leastsqnrm.multiplyenv(pb, ff)

            model_binned = np.zeros((self.fov, self.fov, self.model_over.shape[2]))
            # loop over slices "sl" in the model
            for sl in range(self.model_over.shape[2]):
                model_binned[:, :, sl] = utils.rebin(
                    self.model_over[:, :, sl], (self.over, self.over)
                )

            self.model += w * model_binned

        return self.model

    def fit_image(
        self,
        image,
        reference=None,
        model_in=None,
        savepsfs=False,
        dqm=None,
        weighted=False,
    ):
        """
        Run a least-squares fit on an input image.

        Find the appropriate wavelength scale and rotation.
        If a model is not specified then this
        method will find the appropriate wavelength scale, rotation (and
        hopefully centering as well -- This is not written into the object yet,
        but should be soon).  Without specifying a model, fit_image can take a
        reference image (a cropped deNaNed version of the data) to run
        correlations. It is recommended that the symmetric part of the data be
        used to avoid piston confusion in scaling.

        Parameters
        ----------
        image : 2D float array
            Input image

        reference : 2D float array
            Input reference image

        model_in : 2D array
            Optional model image

        savepsfs : bool
            Save the psfs for writing to file (currently unused)

        dqm : 2D array
            Bad pixel mask of same dimensions as image

        weighted : bool
            Use weighted operations in the least squares routine
        """
        self.model_in = model_in
        self.weighted = weighted
        self.saveval = savepsfs

        if model_in is None:  # No model provided
            # Perform a set of automatic routines
            # A Cleaned up version of your image to enable Fourier fitting for
            # centering crosscorrelation with FindCentering() and
            # magnification and rotation via improve_scaling().

            if reference is None:
                self.reference = image
                if np.isnan(image.any()):
                    raise ValueError(
                        "Must have non-NaN image to "
                        "crosscorrelate for scale. Reference "
                        "image should also be centered."
                    )
            else:
                self.reference = reference

            self.pixscale_measured = self.pixel
            self.fov = image.shape[0]
            self.fittingmodel = self.make_model(
                self.fov,
                bandpass=self.bandpass,
                over=self.over,
                psf_offset=self.bestcenter,
                pixscale=self.pixel,
            )
        else:
            self.fittingmodel = model_in
        if self.weighted is False:
            self.soln, self.residual, self.cond, self.linfit_result = leastsqnrm.matrix_operations(
                image, self.fittingmodel, dqm=dqm
            )
        else:
            self.soln, self.residual, self.cond, self.singvals = leastsqnrm.weighted_operations(
                image, self.fittingmodel, dqm=dqm
            )

        self.rawDC = self.soln[-1]
        self.flux = self.soln[0]
        self.soln = self.soln / self.soln[0]

        # fringephase now in radians
        self.fringeamp, self.fringephase = leastsqnrm.tan2visibilities(self.soln)
        self.fringepistons = utils.fringes2pistons(self.fringephase, len(self.ctrs))
        self.redundant_cps = leastsqnrm.redundant_cps(self.fringephase, n=self.N)
        # RC 8/24
        self.t3_amplitudes = leastsqnrm.t3_amplitudes(self.fringeamp, n=self.N)
        self.redundant_cas = leastsqnrm.closure_amplitudes(self.fringeamp, n=self.N)
        self.q4_phases = leastsqnrm.q4_phases(self.fringephase, n=self.N)  # RC 8/24

    def create_modelpsf(self):
        """Make an image from the object's model and fit solutions by setting modelpsf attribute."""
        try:
            self.modelpsf = np.zeros((self.fov, self.fov))
        except AttributeError:
            self.modelpsf = np.zeros((self.fov_sim, self.fov_sim))

        for ind, coeff in enumerate(self.soln):
            self.modelpsf += self.flux * coeff * self.fittingmodel[:, :, ind]

        return None

    def improve_scaling(self, img):
        """
        Determine the scale and rotation that best fits the data.

        Correlations
        are calculated in the image plane, in anticipation of data with many
        bad pixels.

        Parameters
        ----------
        img : 2D float array
            Input image

        Returns
        -------
        self.pixscale_factor: float
            Improved estimate of pixel scale in radians

        self.rot_measured : float
            Value of mag at the extreme value of rotation from quadratic fit

        self.gof : float
            Goodness of fit
        """
        if not hasattr(self, "bandpass"):
            raise ValueError("This obj has no specified bandpass/wavelength")

        reffov = img.shape[0]
        scal_corrlist = np.zeros((len(self.scallist), reffov, reffov))
        pixscl_corrlist = scal_corrlist.copy()
        scal_corr = np.zeros(len(self.scallist))
        self.pixscl_corr = scal_corr.copy()

        # User can specify a reference set of phases (m) at an earlier point so
        #  that all PSFs are simulated with those phase pistons (e.g. measured
        #  from data at an earlier iteration
        if not hasattr(self, "refphi"):
            self.refphi = np.zeros(len(self.ctrs))
        else:
            pass

        self.pixscales = np.zeros(len(self.scallist))
        for q, scal in enumerate(self.scallist):
            self.test_pixscale = self.pixel * scal
            self.pixscales[q] = self.test_pixscale
            psf = self.simulate(
                bandpass=self.bandpass,
                pixel=self.test_pixscale,
            )
            pixscl_corrlist[q, :, :] = run_data_correlate(img, psf)
            self.pixscl_corr[q] = np.max(pixscl_corrlist[q])
            if True in np.isnan(self.pixscl_corr):
                raise ValueError("Correlation produced NaNs, check your work!")

        self.pixscale_optimal, scal_maxy = utils.findmax(mag=self.pixscales, vals=self.pixscl_corr)
        self.pixscale_factor = self.pixscale_optimal / self.pixel

        radlist = self.rotlist_rad
        corrlist = np.zeros((len(radlist), reffov, reffov))
        self.corrs = np.zeros(len(radlist))

        self.rots = radlist
        for q in range(len(radlist)):
            psf = self.simulate(
                bandpass=self.bandpass,
                fov=reffov,
            )

            corrlist[q, :, :] = run_data_correlate(psf, img)
            self.corrs[q] = np.max(corrlist[q])

        self.rot_measured, maxy = utils.findmax(mag=self.rots, vals=self.corrs)
        self.refpsf = self.simulate(
            bandpass=self.bandpass,
            pixel=self.pixscale_factor * self.pixel,
            fov=reffov,
        )

        try:
            self.gof = goodness_of_fit(img, self.refpsf)
        except Exception:
            self.gof = False

        return self.pixscale_factor, self.rot_measured, self.gof

    def set_pistons(self, phi_m):
        """
        Set piston's phi in meters of OPD at center wavelength LG++.

        Parameters
        ----------
        phi_m : float
            Piston angle
        """
        self.phi = phi_m

    def set_pixelscale(self, pixel_rad):
        """
        Set the detector pixel scale.

        Parameters
        ----------
        pixel_rad : float
            Detector pixel scale
        """
        self.pixel = pixel_rad


def goodness_of_fit(data, bestfit, disk_r=8):
    """
    Calculate goodness of fit between the data and the fit.

    Parameters
    ----------
    data : 2D float array
        Input image

    bestfit : 2D float array
        Fit to input image

    disk_r : int
        Radius of disk

    Returns
    -------
    gof : float
        Goodness of fit
    """
    mask = (
        np.ones(data.shape)
        + utils.makedisk(data.shape[0], 2)
        - utils.makedisk(data.shape[0], disk_r)
    )

    difference = np.ma.masked_invalid(mask * (bestfit - data))

    masked_data = np.ma.masked_invalid(mask * data)

    gof = abs(difference).sum() / abs(masked_data).sum()

    return gof


def run_data_correlate(data, model):
    """
    Calculate correlation between data and model.

    Parameters
    ----------
    data : 2D float array
        Reference image

    model : 2D float array
        Simulated psf

    Returns
    -------
    cor: 2D float array
        Correlation between data and model
    """
    sci = data
    log.debug("shape sci: %s", np.shape(sci))

    cor = utils.rcrosscorrelate(sci, model)

    return cor
