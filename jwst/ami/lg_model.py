# A module for conveniently manipulating an 'NRM object' using the
# Lacour-Greenbaum algorithm. First written by Alexandra Greenbaum in 2014.
import logging
import numpy as np

from . import leastsqnrm as leastsqnrm
from . import analyticnrm2
from . import utils
from . import mask_definitions

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())
log.setLevel(logging.DEBUG)


m = 1.0
mm = 1.0e-3 * m
um = 1.0e-6 * m
mas = 1.0e-3 / (60 * 60 * 180 / np.pi)  # in radians


class NrmModel:
    """
    A class for conveniently dealing with an "NRM object" This should be able
    to take an NRM_mask_definitions object for mask geometry.
    Defines mask geometry and detector-scale parameters.
    Simulates PSF (broadband or monochromatic)
    Builds a fringe model - either by user definition, or automated to data
    Fits model to data by least squares
    Masks: gpi_g10s40, jwst, visir
    Algorithm documented in Greenbaum, A. Z., Pueyo, L. P., Sivaramakrishnan,
    A., and Lacour, S., Astrophysical Journal vol. 798, Jan 2015.
    """

    def __init__(
        self,
        mask=None,
        holeshape="circ",
        pixscale=None,
        over=1,
        pixweight=None,
        datapath="",
        phi=None,
        refdir="",
        chooseholes=False,
        affine2d=None,
        **kwargs,
    ):
        """
        Set attributes of NrmModel class.

        Parameters
        ----------
        mask: string
            keyword for built-in values

        holeshape: string
           shape of apertures

        pixscale: float
           initial estimate of pixel scale in radians

        over: integer
           oversampling factor

        pixweight: 2D float array, default is None
            weighting array

        datapath: string
            directory for output (will remove for final)

        phi: float 1D array
            distance of fringe from hole center in units of waves

        refdir: string
            directory containing ref files (will remove for final)

        chooseholes: list ?
            holes ...?

        affine2d: Affine2d object
            Affine2d object
        """
        if "debug" in kwargs:
            self.debug = kwargs["debug"]
        else:
            self.debug = False

        self.holeshape = holeshape
        self.pixel = pixscale  # det pix in rad (square)
        self.over = over
        self.pixweight = pixweight

        # mask = "jwst"
        # self.maskname = mask

        # get these from mask_definitions instead
        if mask is None:
            log.info("No mask name specified for model, using jwst_g7s6c")
            mask = mask_definitions.NRM_mask_definitions(
                maskname="jwst_g7s6c", chooseholes=chooseholes, holeshape="hex"
            )
        elif isinstance(mask, str):
            mask = mask_definitions.NRM_mask_definitions(
                maskname=mask, chooseholes=chooseholes, holeshape="hex"
            )
        self.ctrs = mask.ctrs
        self.d = mask.hdia
        self.D = mask.activeD

        self.N = len(self.ctrs)
        self.datapath = datapath
        self.refdir = refdir
        self.fmt = "%10.4e"

        # get latest OPD from WebbPSF?

        if phi:  # meters of OPD at central wavelength
            if phi == "perfect":
                self.phi = np.zeros(self.N)  # backwards compatibility
            else:
                self.phi = phi
        else:
            self.phi = np.zeros(self.N)

        self.chooseholes = chooseholes

        # affine2d property not to be changed in NrmModel - create a new
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
        Simulate a detector-scale psf using parameters input from the call and
        already stored in the object,and generate a simulation fits header
        storing all of the  parameters used to generate that psf.  If the input
        bandpass is one number it will calculate a monochromatic psf.

        Parameters
        ----------
        fov: integer, default=None
            number of detector pixels on a side

        bandpass: 2D float array, default=None
            array of the form: [(weight1, wavl1), (weight2, wavl2), ...]

        over: integer
            Oversampling factor

        psf_offset: detector pixels
            Center offset from center of array

        Returns
        -------
        Object's 'psf': float 2D array
            simulated psf
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

    def make_model(
        self, fov=None, bandpass=None, over=1, psf_offset=(0, 0), pixscale=None
    ):
        """
        Generates the fringe model with the attributes of the object using a
        bandpass that is either a single wavelength or a list of tuples of the
        form [(weight1, wavl1), (weight2, wavl2),...].  The model is
        a collection of fringe intensities, where nholes = 7 means the model
        has a @D slice for each of 21 cosines, 21 sines, a DC-like, and a flux
        slice for a toal of 44 2D slices.

        Parameters
        ----------
        fov: integer, default=None
            number of detector pixels on a side

        bandpass: 2D float array, default=None
            array of the form: [(weight1, wavl1), (weight2, wavl2), ...]

        over: integer
           oversampling factor

        psf_offset: detector pixels
            Center offset from center of array

        pixscale: float, default=None
            pixel scale

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

            log.debug("Passed to model_array: psf_offset: {0}".format(psf_offset))
            log.debug("Primary beam in the model created: {0}".format(pb))
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
        pixguess=None,
        rotguess=0,
        psf_offset=(0, 0),
        modelin=None,
        savepsfs=False,
        dqm=None,
        weighted=False,
    ):
        """
        Run a least-squares fit on an input image; find the appropriate
        wavelength scale and rotation. If a model is not specified then this
        method will find the appropriate wavelength scale, rotation (and
        hopefully centering as well -- This is not written into the object yet,
        but should be soon).  Without specifying a model, fit_image can take a
        reference image (a cropped deNaNed version of the data) to run
        correlations. It is recommended that the symmetric part of the data be
        used to avoid piston confusion in scaling.

        Parameters
        ----------
        image: 2D float array
            input image

        reference: 2D float array
            input reference image

        pixguess: float
            estimate of pixel scale of the data

        rotguess: float
            estimate of rotation

        modelin: 2D array
            optional model image

        weighted: boolean
            use weighted operations in the least squares routine

        centering: string, default=None
            type of centering

        savepsfs: boolean
            save the psfs for writing to file (currently unused)

        dqm: 2D array
            bad pixel mask of same dimensions as image

        weighted: boolean
            weight

        Returns
        -------
        None
        """
        self.model_in = modelin
        self.weighted = weighted
        self.saveval = savepsfs

        if modelin is None:  # No model provided
            # Perform a set of automatic routines
            # A Cleaned up version of your image to enable Fourier fitting for
            # centering crosscorrelation with FindCentering() and
            # magnification and rotation via improve_scaling().

            if reference is None:
                self.reference = image
                if np.isnan(image.any()):
                    raise ValueError(
                        "Must have non-NaN image to "
                        + "crosscorrelate for scale. Reference "
                        + "image should also be centered."
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
            self.fittingmodel = modelin
        if self.weighted is False:
            self.soln, self.residual, self.cond, self.linfit_result = (
                leastsqnrm.matrix_operations(image, self.fittingmodel, dqm=dqm)
            )
        else:
            self.soln, self.residual, self.cond, self.singvals = (
                leastsqnrm.weighted_operations(image, self.fittingmodel, dqm=dqm)
            )

        self.rawDC = self.soln[-1]
        self.flux = self.soln[0]
        self.soln = self.soln / self.soln[0]

        # fringephase now in radians
        self.fringeamp, self.fringephase = leastsqnrm.tan2visibilities(self.soln)
        self.fringepistons = utils.fringes2pistons(self.fringephase, len(self.ctrs))
        self.redundant_cps = leastsqnrm.redundant_cps(self.fringephase, n=self.N)
        self.redundant_cas = leastsqnrm.closure_amplitudes(self.fringeamp, n=self.N)

    def create_modelpsf(self):
        """
        Make an image from the object's model and fit solutions, by setting the
        NrmModel object's modelpsf attribute

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        try:
            self.modelpsf = np.zeros((self.fov, self.fov))
        except AttributeError:
            self.modelpsf = np.zeros((self.fov_sim, self.fov_sim))

        for ind, coeff in enumerate(self.soln):
            self.modelpsf += self.flux * coeff * self.fittingmodel[:, :, ind]

        return None

    def improve_scaling(
        self, img, scaleguess=None, rotstart=0.0, centering="PIXELCENTERED"
    ):
        """
        Determine the scale and rotation that best fits the data.  Correlations
        are calculated in the image plane, in anticipation of data with many
        bad pixels.

        Parameters
        ----------
        img: 2D float array
            input image

        scaleguess: float
            initial estimate of pixel scale in radians

        rotstart: float
            estimate of rotation

        centering: string, default='PIXELCENTERED'
            type of centering

        Returns
        -------
        self.pixscale_factor: float
            improved estimate of pixel scale in radians

        self.rot_measured: float
            value of mag at the extreme value of rotation from quadratic fit

        self.gof: float
            goodness of fit
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

        self.pixscale_optimal, scal_maxy = utils.findmax(
            mag=self.pixscales, vals=self.pixscl_corr
        )
        self.pixscale_factor = self.pixscale_optimal / self.pixel

        radlist = self.rotlist_rad
        corrlist = np.zeros((len(radlist), reffov, reffov))
        self.corrs = np.zeros(len(radlist))

        self.rots = radlist
        for q, rad in enumerate(radlist):
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
        Set piston's phi in meters of OPD at center wavelength LG++

        Parameters
        ----------
        phi: float
            piston angle

        Returns
        -------
        None

        """
        self.phi = phi_m

    def set_pixelscale(self, pixel_rad):
        """
        Set the detector pixel scale

        Parameters
        ----------
        pixel_rad: float
            Detector pixel scale

        Returns
        -------
        None
        """
        self.pixel = pixel_rad


def goodness_of_fit(data, bestfit, diskR=8):
    """
    Calculate goodness of fit between the data and the fit.

    Parameters
    ----------
    data: 2D float array
        input image

    bestfit: 2D float array
        fit to input image

    diskR: integer
        radius of disk

    Returns
    -------
    gof: float
        goodness of fit
    """
    mask = (
        np.ones(data.shape)
        + utils.makedisk(data.shape[0], 2)
        - utils.makedisk(data.shape[0], diskR)
    )

    difference = np.ma.masked_invalid(mask * (bestfit - data))

    masked_data = np.ma.masked_invalid(mask * data)

    gof = abs(difference).sum() / abs(masked_data).sum()

    return gof


def run_data_correlate(data, model):
    """
    Calculate correlation between data and model

    Parameters
    ----------
    data: 2D float array
        reference image

    model: 2D float array
        simulated psf

    Returns
    -------
    cor: 2D float array
        correlation between data and model

    """
    sci = data
    log.debug("shape sci: %s", np.shape(sci))

    cor = utils.rcrosscorrelate(sci, model)

    return cor
