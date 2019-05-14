#
# A module for conveniently manipulating an 'NRM object' using the
# Lacour-Greenbaum algorithm. First written by Alexandra Greenbaum in 2014.
#
# This module:
#   Defines mask geometry and detector-scale parameters
#   Simulates PSF (broadband or monochromatic)
#   Builds a fringe model - either by user definition, or automated to data
#   Fits model to data by least squares
#
# Algorithm documented in:  Greenbaum, A. Z., Pueyo, L. P.,
# Sivaramakrishnan, A., and Lacour, S. ; Astrophysical Journal (submitted) 2014.
# Developed with NASA APRA (AS, AZG), NSF GRFP (AZG), NASA Sagan (LP), and
# French taxpayer (SL) support.
#
# Heritage mathematica nb from Alex Greenbaum & Laurent Pueyo
# Heritage python by Alex Greenbaum & Anand Sivaramakrishnan Jan 2013
# - updated May 2013 to include hexagonal envelope
# - updated (hard refactored) Oct-Nov 2014 Anand S.

import logging
import numpy as np

from . import leastsqnrm as leastsqnrm
from . import analyticnrm2
from . import utils
from . import hexee
from . import nrm_consts

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class NrmModel:

    def __init__(self, mask=None, holeshape="circ", pixscale=hexee.mas2rad(65),
                 rotate=False, over=1, flip=False, pixweight=None, scallist=None,
                 rotlist_deg=None, phi="perfect"):
        """
        Short Summary
        -------------
        Set attributes of NrmModel class.

        Parameters
        ----------
        mask: string
            keyword for built-in values

        holeshape: string
           shape of apertures

        pixscale: float
           initial estimate of pixel scale in radians

        rotate: float
           initial estimate of rotation in radians

        over: integer
           oversampling factor

        flip: Boolean, default=False
            change sign of 2nd coordinate of holes

        pixweight: 2D float array, default is None
            weighting array

        scallist: float 1D array
           candidate relative pixel scales

        rotlist_deg: float 1D array
            Search window for rotation fine-tuning, in degrees

        phi: float 1D array
            distance of fringe from hole center in units of waves
        """
        self.holeshape = holeshape
        self.pixel = pixscale
        self.over = over
        self.maskname = mask
        self.pixweight = pixweight

        if mask.lower() == 'jwst':
            self.ctrs = np.array( [[ 0.00000000, -2.640000],
                                   [-2.2863100, 0.0000000],
                                   [ 2.2863100 , -1.3200001],
                                   [-2.2863100 ,  1.3200001],
                                   [-1.1431500 ,  1.9800000],
                                   [ 2.2863100 ,  1.3200001],
                                   [ 1.1431500 ,  1.9800000]] )
            self.d = 0.80
            self.D = 6.5
        else:
            try:
                log.debug('mask.ctrs:%s', mask.ctrs)
            except AttributeError:
                raise AttributeError("mask must be either 'jwst' \
                                                      or NRM_mask_geometry object")

        log.debug('NrmModel: ctrs flipped in init for CV1, CV2')

        if rotate:
            log.info('Providing additional rotation %s degrees',
                     rotate * 180. / np.pi)

            # Rotates vector counterclockwise in coordinates
            self.rotation = rotate
            self.ctrs = leastsqnrm.rotatevectors(self.ctrs, rotate)

            # From now on this 'rotated' set of centers is used as the
            #  nominal, and rotation searches (using rotlist_rad) are
            #  performed with this rotated version of the 'as designed'
            #  mask..  In CV1 and CV2 the mask is "flipped" by
            #  multiplying ctrs[:1] by -1... which places segment B4
            #  (the 6 o clock segment) at the top instead of at bottom
            #  in traditional XY plots

        self.N = len(self.ctrs)

        if scallist is None:
            self.scallist = np.array([0.995, 0.998, 1.0, 1.002, 1.005, ])
        else:
            self.scallist = scallist

        if rotlist_deg is None:
            self.rotlist_rad = np.array([-1.0,-0.5,0.0,0.5,1.0]) * np.pi / 180.0
        else:
            self.rotlist_rad = rotlist_deg * np.pi / 180.0

        if phi == "perfect":
            self.phi = np.zeros(len(self.ctrs))
        elif phi == 'nb':
            self.phi = nrm_consts.phi_nb
        else:
            self.phi = phi


    def simulate(self, fov=None, bandpass=None, over=None, pixweight=None,
                 pixel=None, rotate=False, centering="PIXELCENTERED"):
        """
        Short Summary
        -------------
        Simulate a psf using parameters input from the call and already stored in
        the object. It also generates a simulation fits header storing all of the
        parameters used to generate that psf.  If the input bandpass is one
        number it will calculate a monochromatic PSF.

        Parameters
        ----------
        fov: integer, default=None
            number of detector pixels on a side

        bandpass: 2D float array, default=None
            array of the form: [(weight1, wavl1), (weight2, wavl2), ...]

        over: integer
            Oversampling factor

        pixweight: 2D float array, default=None
            weighting array

        pixel: float, default=None
            pixel scale

        rotate: float, default=False,
            rotation angle in radians

        centering: string, default=None
            type of centerings

        Returns
        -------
        Object's 'psf': float 2D array
            simulated psf
        """
        # First set up conditions for choosing various parameters
        if fov is None:
            if not hasattr(self, 'fov'):
                log.critical('Field is not specified')
                return None
            else:
                self.fov_sim = self.fov
                log.debug('Using predefined FOV size: %s', self.fov)
        else:
            self.fov_sim = fov

        if hasattr(centering, '__iter__'):
            if centering == 'PIXELCENTERED':
                centering=(0.5, 0.5)
            elif centering == 'PIXELCORNER':
                centering=(0.0, 0.0)

        self.bandpass = bandpass

        if not hasattr(self, 'over'):
            if over is None:
                self.over = 1
            else:
                self.over = over

        if self.over is None:
            self.over = over

        if pixweight is not None:
            self.over = self.pixweight.shape[0]

        self.phi = np.zeros(len(self.ctrs))

        if rotate: # this is a 'temporary' rotation of self.ctrs
            #  without altering self.ctrs
            self.rotate = rotate
            self.rotctrs = leastsqnrm.rotatevectors(self.ctrs, self.rotate)
        else:
            self.rotctrs = self.ctrs

        if pixel is None:
            self.pixel_sim = self.pixel
        else:
            self.pixel_sim = pixel

        # The polychromatic case:
        if hasattr(self.bandpass, '__iter__'):
            log.debug("------Simulating Polychromatic------")
            self.psf_over = np.zeros((self.over*self.fov_sim,
                                      self.over*self.fov_sim))
            for w,l in self.bandpass: # w: weight, l: lambda (wavelength)
                self.psf_over += w*analyticnrm2.PSF(self.pixel_sim,
                            self.fov_sim, self.over, self.rotctrs, self.d, l,
                            self.phi, centering = centering, shape=self.holeshape)

            log.debug("BINNING UP TO PIXEL SCALE")

        # The monochromatic case if bandpass input is a single wavelength
        else:
            self.lam = bandpass

            log.debug("Calculating Oversampled PSF")
            self.psf_over = analyticnrm2.PSF(self.pixel_sim, self.fov_sim,
                            self.over, self.rotctrs, self.d, self.lam,
                            self.phi, centering=centering,
                            shape=self.holeshape)

        self.psf = utils.rebin(self.psf_over, (self.over, self.over))

        return self.psf


    def make_model(self, fov=None, bandpass=None, over=False,
                    centering='PIXELCENTERED', pixweight=None, pixscale=None,
                    rotate=False, flip=False):
        """
        Short Summary
        -------------
        Generates the fringe model with the attributes of the object
        using a bandpass as a list of tuples.

        Parameters
        ----------
        fov: integer, default=None
            number of detector pixels on a side

        bandpass: 2D float array, default=None
            array of the form: [(weight1, wavl1), (weight2, wavl2), ...]

        over: integer
           oversampling factor

        centering: string, default=None
            type of centering

        pixweight: 2D float array, default=None
            weighting array

        pixscale: float, default=None
            pixel scale

        rotate: float, default=False
            rotation angle in radians

        flip: Boolean, default=False
            change sign of 2nd coordinate of holes

        Returns
        -------
        Object's 'model': fringe model
            Generated fringe model
        """
        if fov:
            self.fov = fov

        if over is False:
            self.over = 1
        else:
            self.over = over

        if pixweight is not None:
            self.over = self.pixweight.shape[0]

        if hasattr(self, 'pixscale_measured'):
            if self.pixscale_measured is not None:
                self.modelpix = self.pixscale_measured

        if pixscale is None:
            self.modelpix = self.pixel
        else:
            self.modelpix = pixscale

        if rotate:
            if flip is True:
                self.modelctrs = leastsqnrm.flip(
                    leastsqnrm.rotatevectors(self.ctrs, self.rot_measured))
            else:
                self.modelctrs = leastsqnrm.rotatevectors(
                    self.ctrs, self.rot_measured)
        else:
            self.modelctrs = self.ctrs

        if not hasattr(bandpass, '__iter__'):
            self.lam = bandpass
            self.model = np.ones((self.fov, self.fov, self.N*(self.N-1)+2))
            self.model_beam, self.fringes = leastsqnrm.model_array(
                     self.modelctrs, self.lam, self.over, self.modelpix,
                     self.fov, self.d, shape=self.holeshape, centering=centering)

            log.debug("centering: {0}".format(centering))
            log.debug("what primary beam has the model created?"+
                                " {0}".format(self.model_beam))

            # this routine multiplies the envelope by each fringe "image"
            self.model_over = leastsqnrm.multiplyenv(self.model_beam, self.fringes)

            self.model = np.zeros((self.fov,self.fov, self.model_over.shape[2]))
            # loop over slices "sl" in the model
            for sl in range(self.model_over.shape[2]):
                self.model[:,:,sl] = utils.rebin( self.model_over[:,:,sl],
                                                (self.over, self.over))

            return self.model

        else:
            self.bandpass = bandpass

            # The model shape is (fov) x (fov) x (# solution coefficients)
            # the coefficient refers to the terms in the analytic equation
            # There are N(N-1) independent pistons, double-counted by cosine
            # and sine, one constant term and a DC offset.
            self.model = np.ones((self.fov, self.fov, self.N*(self.N-1)+2))
            self.model_beam = np.zeros((self.over*self.fov, self.over*self.fov))
            self.fringes = np.zeros((
                self.N*(self.N-1)+1, self.over*self.fov, self.over*self.fov))

            for w,l in self.bandpass: # w: weight, l: lambda (wavelength)
                # model_array returns the envelope and fringe model
                pb, ff = leastsqnrm.model_array(
                    self.modelctrs, l, self.over, self.modelpix, self.fov,
                    self.d, shape=self.holeshape, centering=centering)

                log.debug("centering: {0}".format(centering))
                log.debug("what primary beam has the model created? {0}".format(pb))

                self.model_beam += pb
                self.fringes += ff

                # this routine multiplies the envelope by each fringe "image"
                self.model_over = leastsqnrm.multiplyenv(pb, ff)

                model_binned = np.zeros((
                    self.fov,self.fov, self.model_over.shape[2]))

                # loop over slices "sl" in the model
                for sl in range(self.model_over.shape[2]):
                    model_binned[:,:,sl] = utils.rebin(
                        self.model_over[:,:,sl], (self.over, self.over))

                    self.model += w*model_binned

            return self.model


    def fit_image(self, image, reference=None, pixguess=None, rotguess=0,
                  modelin=None, weighted=False, centering='PIXELCENTERED',
                  savepsfs=True):
        """
        Short Summary
        -------------
        Run a least-squares fit on an input image; find the appropriate
        wavelength scale and rotation. If a model is not specified then this
        method will find the appropriate wavelength scale, rotation (and
        hopefully centering as well -- This is not written into the object yet,
        but should be soon).

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

        Returns
        -------
        None
        """
        self.model_in = modelin
        self.weighted = weighted
        self.saveval = savepsfs

        if modelin is None: # No model provided
            # Perform a set of automatic routines
            # A Cleaned up version of your image to enable Fourier fitting for
            # centering crosscorrelation with FindCentering() and
            # magnification and rotation via improve_scaling().

            if reference is None:
                self.reference = image
                if np.isnan(image.any()):
                    raise ValueError("Must have non-NaN image to "+
                        "crosscorrelate for scale. Reference "+
                        "image should also be centered. Get to it.")
            else:
                self.reference = reference

            if pixguess is None or rotguess is None:
                raise ValueError("MUST SPECIFY GUESSES FOR PIX & ROT")

            self.improve_scaling(self.reference, scaleguess=self.pixel,
                rotstart=rotguess, centering=centering)

            self.pixscale_measured = self.pixscale_factor*self.pixel

            self.fov = image.shape[0]
            self.fittingmodel = self.make_model(self.fov, bandpass=self.bandpass,
                over=self.over, rotate=True, centering=centering,
                pixscale=self.pixscale_measured)
        else:
            self.fittingmodel = modelin

        if weighted is not False:
            self.soln, self.residual = leastsqnrm.weighted_operations(image,
                self.fittingmodel, weights=self.weighted)
        else:
            self.soln, self.residual, self.cond = leastsqnrm.matrix_operations(
                image, self.fittingmodel)

        self.rawDC = self.soln[-1]
        self.flux = self.soln[0]
        self.soln = self.soln/self.soln[0]
        self.deltapsin = leastsqnrm.sin2deltapistons(self.soln)
        self.deltapcos = leastsqnrm.cos2deltapistons(self.soln)

        self.fringeamp, self.fringephase = leastsqnrm.tan2visibilities(self.soln)
        self.piston = utils.fringes2pistons(self.fringephase, len(self.ctrs))
        self.closurephase = leastsqnrm.closurephase(self.fringephase, N=self.N)
        self.redundant_cps = leastsqnrm.redundant_cps(self.fringephase, N=self.N)
        self.redundant_cas = leastsqnrm.return_CAs(self.fringeamp, N=self.N)


    def create_modelpsf(self):
        """
        Short Summary
        -------------
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


    def improve_scaling(self, img, scaleguess=None, rotstart=0.0,
                        centering='PIXELCENTERED'):
        """
        Short Summary
        -------------
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
        if not hasattr(self, 'bandpass'):
            raise ValueError("This object has no specified bandpass/wavelength")

        reffov = img.shape[0]
        scal_corrlist = np.zeros((len(self.scallist), reffov, reffov))
        pixscl_corrlist = scal_corrlist.copy()
        scal_corr = np.zeros(len(self.scallist))
        self.pixscl_corr = scal_corr.copy()

        # User can specify a reference set of phases (m) at an earlier point so
        #  that all PSFs are simulated with those phase pistons (e.g. measured
        #  from data at an earlier iteration
        if not hasattr(self, 'refphi'):
            self.refphi = np.zeros(len(self.ctrs))
        else:
            pass

        self.pixscales = np.zeros(len(self.scallist))
        for q, scal in enumerate(self.scallist):
            self.test_pixscale = self.pixel*scal
            self.pixscales[q] = self.test_pixscale
            psf = self.simulate(bandpass=self.bandpass, fov=reffov,
                    pixel = self.test_pixscale, centering=centering)
            pixscl_corrlist[q,:,:] = run_data_correlate(img,psf)
            self.pixscl_corr[q] = np.max(pixscl_corrlist[q])
            if True in np.isnan(self.pixscl_corr):
                raise ValueError("Correlation produced NaNs, check your work!")

        self.pixscale_optimal, scal_maxy = utils.findmax(
            mag=self.pixscales, vals=self.pixscl_corr)
        self.pixscale_factor = self.pixscale_optimal / self.pixel

        radlist = self.rotlist_rad
        corrlist = np.zeros((len(radlist), reffov, reffov))
        self.corrs = np.zeros(len(radlist))

        self.rots = radlist
        for q,rad in enumerate(radlist):
            psf = self.simulate(bandpass=self.bandpass, fov=reffov,
                    pixel=self.pixscale_optimal, rotate=rad, centering=centering)

            corrlist[q,:,:] = run_data_correlate(psf, img)
            self.corrs[q] = np.max(corrlist[q])

        self.rot_measured, maxy = utils.findmax(mag=self.rots, vals = self.corrs)
        self.refpsf = self.simulate(bandpass=self.bandpass,
                pixel=self.pixscale_factor*self.pixel, fov=reffov,
                rotate=self.rot_measured, centering=centering)

        try:
            self.gof = goodness_of_fit(img, self.refpsf)
        except Exception:
            self.gof = False

        return self.pixscale_factor, self.rot_measured, self.gof


def makedisk(N, R, ctr=(0,0)):
    """
    Short Summary
    -------------
    Calculate a 'disk', an array whose values =1 in a circular region near
    the center of the array, and =0 elsewhere. (Anand's emailed version)

    Parameters
    ----------
    N: integer
        size of 1 dimension of the array to be returned

    R: integer
        radius of disk

    ctr: (integer, integer)
        center of disk

    array: 'ODD' or 'EVEN'
        parity of size of edge

    Returns
    -------
    array: 2D integer array
    """
    if N%2 == 1: # odd
        M = (N-1)/2
        xx = np.linspace(-M-ctr[0],M-ctr[0],N)
        yy = np.linspace(-M-ctr[1],M-ctr[1],N)
    if N%2 == 0: # even
        M = N/2
        xx = np.linspace(-M-ctr[0],M-ctr[0]-1,N)
        yy = np.linspace(-M-ctr[1],M-ctr[1]-1,N)

    (x,y) = np.meshgrid(xx, yy.T)
    r = np.sqrt((x**2)+(y**2))
    array = np.zeros((N,N))
    array[r<R] = 1

    return array


def goodness_of_fit(data, bestfit, diskR=8):
    """
    Short Summary
    -------------
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
    mask = np.ones(data.shape) + makedisk(data.shape[0], 2) -\
                    makedisk(data.shape[0], diskR)

    difference = np.ma.masked_invalid(mask * (bestfit - data))

    masked_data = np.ma.masked_invalid(mask * data)

    gof = abs(difference).sum() / abs(masked_data).sum()

    return gof


def run_data_correlate(data, model):
    """
    Short Summary
    -------------
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
    log.debug('shape sci: %s', np.shape(sci))

    cor = utils.rcrosscorrelate(sci, model)

    return cor
