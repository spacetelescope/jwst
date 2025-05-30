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
    A class for conveniently dealing with an NRM object.

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
        pixscale,
        bandpass,
        mask="jwst_ami",
        holeshape="hex",
        over=1,
        phi=None,
        chooseholes=None,
        affine2d=None,
    ):
        """
        Set attributes of LgModel class.

        Parameters
        ----------
        nrm_model : NRMModel (or NRMDefinition object? test.)
            Datamodel containing mask geometry information
        pixscale : float
           Initial estimate of pixel scale in radians
        bandpass : np.ndarray[float]
            Array of the form: [(weight1, wavl1), (weight2, wavl2), ...]
        mask : str
            Keyword for built-in values
        holeshape : str
           Shape of apertures, default="hex"
        over : int
           Oversampling factor
        phi : float 1D array
            Distance of fringe from hole center in units of waves
        chooseholes : list of strings, default None
            E.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask
            If None, use the real seven-hole mask.
        affine2d : Affine2d object
            Affine2d object
        """
        self.bandpass = bandpass
        self.holeshape = holeshape
        self.pixel = pixscale  # det pix in rad (square)
        self.over = over

        self.mask = mask_definition_ami.NRMDefinition(
            nrm_model, maskname=mask, chooseholes=chooseholes
        )
        self.ctrs = self.mask.ctrs
        self.d = self.mask.hdia
        self.D = self.mask.active_D

        self.N = len(self.ctrs)
        self.fmt = "%10.4e"

        # get closest in time OPD from STPSF?

        if (phi is None) or (phi == "perfect"):  # meters of OPD at central wavelength
            self.phi = np.zeros(self.N)  # backwards compatibility
        else:
            self.phi = phi

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

    def simulate(self, fov, psf_offset=(0, 0)):
        """
        Simulate a detector-scale psf.

        Use parameters input from the call and
        already stored in the object, and generate a simulation fits header
        storing all of the  parameters used to generate that psf.  If the input
        bandpass is one number it will calculate a monochromatic psf.

        Parameters
        ----------
        fov : int
            Number of detector pixels on a side
        psf_offset : detector pixels
            Center offset from center of array

        Returns
        -------
        Object's 'psf' : float 2D array
            Simulated psf of shape (fov, fov).
        """
        # First set up conditions for choosing various parameters

        self.psf_over = np.zeros((self.over * fov, self.over * fov))
        nspec = 0
        # accumulate polychromatic oversampled psf in the object

        for w, l in self.bandpass:  # w: wavelength's weight, l: lambda (wavelength)
            self.psf_over += w * analyticnrm2.psf(
                self.pixel,  # det pixel, rad
                fov,  # in detpix number
                self.over,
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
        self.psf = utils.rebin(self.psf_over, (self.over, self.over))

        return self.psf

    def make_model(self, fov, psf_offset=(0, 0)):
        """
        Generate the fringe model.

        Use the attributes of the object with a bandpass that is either a single
        wavelength or a list of tuples of the form
        [(weight1, wavl1), (weight2, wavl2),...].  The model is
        a collection of fringe intensities, where nholes = 7 means the model
        has a @D slice for each of 21 cosines, 21 sines, a DC-like, and a flux
        slice for a total of 44 2D slices.

        Parameters
        ----------
        fov : int
            Number of detector pixels on a side
        psf_offset : detector pixels
            Center offset from center of array

        Returns
        -------
        model : np.ndarray[float]
            Generated fringe model, shape (fov, fov, N * (N - 1) + 2)
        """
        self.fov = fov

        # The model shape is (fov) x (fov) x (# solution coefficients)
        # the coefficient refers to the terms in the analytic equation
        # There are N(N-1) independent pistons, double-counted by cosine
        # and sine, one constant term and a DC offset.

        self.model = np.zeros((self.fov, self.fov, self.N * (self.N - 1) + 2))
        self.model_beam = np.zeros((self.over * self.fov, self.over * self.fov))
        self.fringes = np.zeros(
            (self.N * (self.N - 1) + 1, self.over * self.fov, self.over * self.fov)
        )

        for w, l in self.bandpass:  # w: weight, l: lambda (wavelength)
            # model_array returns the envelope and fringe model as a list of
            #   oversampled fov x fov slices
            pb, ff = analyticnrm2.model_array(
                self.ctrs,
                l,
                self.over,
                self.pixel,
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
            model_over = leastsqnrm.multiplyenv(pb, ff)

            model_binned = np.zeros((self.fov, self.fov, model_over.shape[2]))
            # loop over slices "sl" in the model
            for sl in range(model_over.shape[2]):
                model_binned[:, :, sl] = utils.rebin(model_over[:, :, sl], (self.over, self.over))

            self.model += w * model_binned

        return self.model

    def fit_image(
        self,
        image,
        model_in,
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
        TODO: change name of self.singvals or self.linfit_results to be the same, for consistency.
        This would be easier if matrix_operations and weighted_operations both did their fitting
        with scipy

        Parameters
        ----------
        image : 2D float array
            Input image
        model_in : 2D float array
            Model image
        dqm : 2D array
            Bad pixel mask of same dimensions as image
        weighted : bool
            Use weighted operations in the least squares routine
        """
        self.weighted = weighted
        self.fittingmodel = model_in
        if dqm is None:
            dqm = np.zeros(image.shape, dtype="bool")

        if not weighted:
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
    float
        Goodness of fit
    """
    mask = (
        np.ones(data.shape)
        + utils.makedisk(data.shape[0], 2)
        - utils.makedisk(data.shape[0], disk_r)
    )

    difference = np.ma.masked_invalid(mask * (bestfit - data))

    masked_data = np.ma.masked_invalid(mask * data)

    return abs(difference).sum() / abs(masked_data).sum()


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
    2D float array
        Correlation between data and model
    """
    sci = data
    log.debug("shape sci: %s", np.shape(sci))

    return utils.rcrosscorrelate(sci, model)
