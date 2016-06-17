#
# A module for conveniently manipulating an 'NRM object' using the
# Lacour-Greenbaum algorithm. First written by Alexandra Greenbaum in 2014.
# Algorithm documented in:  Greenbaum, A. Z., Pueyo, L. P.,
# Sivaramakrishnan, A., and Lacour, S. ; Astrophysical Journal (submitted) 2014.
# Developed with NASA APRA (AS, AZG), NSF GRFP (AZG), NASA Sagan (LP), and
# French taxpayer (SL) support.
#
# Heritage mathematica nb from Alex Greenbaum & Laurent Pueyo
# Heritage python by Alex Greenbaum & Anand Sivaramakrishnan Jan 2013
# - updated May 2013 to include hexagonal envelope
# - updated (hard refactored) Oct-Nov 2014 Anand S.

from __future__ import absolute_import, division

import logging
import numpy as np
from . import leastsqnrm as leastsqnrm
from . import analyticnrm2
from . import utils
from . import NRM_consts

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class NRM_Model():

    def __init__(self, mask=None, holeshape="circ",
                  pixscale=leastsqnrm.mas2rad(64), rotate=False,
                  over=None, scallist=None, rotlist_deg=None):
        """
        Short Summary
        -------------
        Set attributes of NRM_Model object.

        Parameters
        ----------
        mask: string
            keyword for built-in values

        holeshape: string
           shape of apertures

        pixscale: float
           initial estimate of pixel scale in radians.??

        rotate: float
           initial estimate of rotation in radians.??

        over: integer
           oversampling factor

        scallist: float 1D array
           candidate relative pixel scales

        rotlist_deg: float 1D array
            Search window for rotation fine-tuning, in degrees

        """
        self.holeshape = holeshape
        self.pixel = pixscale

        if over:
            self.over = over
        else:
            self.over = 1

        self.maskname = mask
        self.ctrs = NRM_consts.ctrs # this is in the 'flip' code block
        self.ctrs = leastsqnrm.flip(self.ctrs)

        self.d = 0.80  # dg - what is this ?
        self.D = 6.5   # dg - what is this ?
        self.N = len(self.ctrs)  # dg - what is this ?

        log.info('NRM_Model: ctrs flipped in init for CV1, CV2')

        if rotate: # incoming, in radians
            log.info('Providing additional rotation %s degrees',
                     rotate * 180. / np.pi)

            # Rotates vector counterclockwise in coordinates
            self.ctrs = leastsqnrm.rotatevectors(self.ctrs, rotate)

            # from now on this 'rotated' set of centers is used as the
            #  nominal, and rotation searches (using rotlist_rad) are
            #  performed with this rotated version of the 'as designed'
            #  mask..  In CV1 and CV2 the mask is "flipped" by
            # multiplying ctrs[:1] by -1... which places segment B4
            # (the 6 o clock segment)
            # at the top instead of at bottom in traditional XY plots

        if scallist is None:
            self.scallist = np.array([0.995, 0.998, 1.0, 1.002, 1.005, ])
        else:
            self.scallist = scallist

        if rotlist_deg is None:
            self.rotlist_rad = np.array([0.0, 0.5, 1.0, 1.5, 2.0]) *\
                               np.pi / 180.0
        else:
            self.rotlist_rad = rotlist_deg * np.pi / 180.0


    def simulate(self, fov=None, bandpass=None, pixel=None, rotate=False,
                  centering=None):
        """
        Short Summary
        -------------
        Simulate a psf using parameters input from the call and already
        stored in the object.

        Parameters
        ----------
        fov: integer, default=None
            number of detector pixels on a side

        bandpass: 2D float array, default=None
            array of the form: [(weight1, wavl1), (weight2, wavl2), ...]

        pixel: float, default=None
            pixel scale

        rotate: float, default=False,
            rotation angle in radians

        centering: string, default=None
            type of centering

        Returns
        -------
        Object's 'psf': float 2D array
            simulated psf
        """
        # First set up conditions for choosing various parameters
        if fov == None:
            if not hasattr(self, 'fov'):
                log.critical('You must specify a field')
            else:
                self.fov_sim = self.fov
                log.debug('Using predefined FOV size: %s', self.fov)
        else:
            self.fov_sim = fov

        if centering is not None:
            self.centering = centering

        if not hasattr(self, 'centering'):  # dg - perhaps set this earlier
            self.centering = 'PIXELCENTERED'

        self.bandpass = bandpass  #  dg - perhaps set this in __init__

        if not hasattr(self, 'over'):
            if over == None:
                self.over = 1
            else:
                self.over = over

        if self.over == None:
            self.over = over

        self.phi = np.zeros(len(self.ctrs))

        if rotate: # this is a 'temporary' rotation of self.ctrs
            #  without altering self.ctrs
            self.rotate = rotate
            self.rotctrs = leastsqnrm.rotatevectors(self.ctrs, self.rotate)
        else:
            self.rotctrs = self.ctrs # if initialization included a
            #   rotation self.ctrs has been rotated

        if pixel == None:
            self.pixel_sim = self.pixel
        else:
            self.pixel_sim = pixel

        # The polychromatic case:
        log.debug("------Simulating Polychromatic------")

        self.psf_over = np.zeros((self.over * self.fov_sim,
                        self.over * self.fov_sim))

        for w, l in self.bandpass: # w: weight, l: lambda (wavelength)
            log.debug('weight: %s', w)
            log.debug('wavelength: %s', l)
            log.debug('fov: %s', self.fov_sim)
            log.debug('over: %s', self.over)

            self.psf_over += w * analyticnrm2.PSF(self.pixel_sim, \
                        self.fov_sim, self.over, self.rotctrs, self.d, l, \
                        self.phi, centering=self.centering, \
                        shape=self.holeshape)

        log.debug("BINNING UP TO PIXEL SCALE")

        self.psf = utils.rebin(self.psf_over, (self.over, self.over))
        self.psf = self.psf[::-1, :]

        return self.psf


    def make_model(self, fov=None, bandpass=None, over=False,
                    centering='PIXELCENTERED', pixscale=None, rotate=False):
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

        pixscale: float, default=None
            pixel scale

        rotate: float, default=False
            rotation angle in radians

        Returns
        -------
        Object's 'model': fringe model
            Generated fringe model
        """

        if self.maskname.lower() == 'jwst':
            self.ctrs = leastsqnrm.flip(self.ctrs)

        if fov:
            self.fov = fov

        if over is False:
            self.over = 1
        else:
            self.over = over

        if not hasattr(self, 'centering'):
            self.centering = centering

        if hasattr(self, 'pixscale_measured'):
            if self.pixscale_measured is not None:
                self.modelpix = self.pixscale_measured

        if pixscale == None:
            self.modelpix = self.pixel
        else:
            self.modelpix = pixscale

        if rotate:
            log.debug('ROTATE: %s', rotate)
            self.modelctrs = leastsqnrm.rotatevectors(self.ctrs,
                                                      self.rot_measured)
            log.debug('MODEL USING MEASURED ROTATION')

        log.debug('So the bandpass should have shape: %s', np.shape(bandpass))
        log.debug('And should be iterable: %s', hasattr(bandpass, '__iter__'))

        self.bandpass = bandpass

        log.debug('-------Polychromatic bandpass---------')
        log.debug('%s', self.bandpass)

        # The model shape is (fov) x (fov) x (# solution coefficients)
        # the coefficient refers to the terms in the analytic equation.
        # There are N(N-1) independent pistons, double-counted by cosine
        # and sine, one constant term and a DC offset.
        self.model = np.ones((self.fov, self.fov, self.N * (self.N - 1) + 2))
        self.model_beam = np.zeros((self.over * self.fov, self.over * self.fov))
        self.fringes = np.zeros((self.N * (self.N - 1) + 1, \
                    self.over * self.fov, self.over * self.fov))
        for w, l in self.bandpass: # w: weight, l: lambda (wavelength)
            log.debug('weight:%s lambda: %s', w, l)

            # model_array returns the envelope and fringe model
            pb, ff = leastsqnrm.model_array(self.modelctrs, l, self.over, \
                              self.modelpix, self.fov, self.d, \
                              shape=self.holeshape, centering=self.centering)

            log.debug("centering: {0}".format(self.centering))
            log.debug("what primary beam has model" + "created? {0}".format(pb))

            self.model_beam += pb
            self.fringes += ff

            # this routine multiplies the envelope by each fringe "image"
            self.model_over = leastsqnrm.multiplyenv(pb, ff)

            log.debug('NRM MODEL model shape: %s', self.model_over.shape)

            model_binned = np.zeros((self.fov, self.fov,
                        self.model_over.shape[2]))

            # loop over slices "sl" in the model
            for sl in range(self.model_over.shape[2]):
                model_binned[:, :, sl] = utils.rebin(self.model_over[:, :, sl], \
                (self.over, self.over))

            self.model += w * model_binned

        return self.model


    # dg - leave in 'modelin' and 'centering' for now
    def fit_image(self, image, pixguess=None, modelin=None,
                   centering='PIXELCENTERED'):
        """
        Short Summary
        -------------
        Run a least-squares fit on an input image.
        Specifying a model is optional. If a model is not specified then this
        method will find the appropriate wavelength scale, rotation (and
        hopefully centering as well -- This is not written into the object
        yet, but should be soon).

        Without specifying a model, fit_image can take a reference image
        (a cropped deNaNed version of the data) to run correlations. It is
        recommended that the symmetric part of the data be used to avoid piston
        confusion in scaling.

        Parameters
        ----------
        image: 2D float array
            input image

        pixguess: float
            estimate of pixel scale of the data

        modelin: 2D array
            optional model image

        centering: string, default=None
             type of centering

        Returns
        -------
        None

        """
        self.model_in = modelin

        if hasattr(self, 'lam'):
            self.bandpass = self.lam

        if modelin is None:
            # dg - find out from Anand if I can delete this block

            if pixguess == None:
                log.critical('YOU MUST SPECIFY GUESSES FOR PIX & ROT.')
                log.critical('fit image requires a pixel scale guess PIXSCALE')
                log.critical('  and rotation guess ROTGUESS keyword')

            # A Cleaned up version of your image to enable Fourier fitting for
            # magnification and rotation via refine_scaling().

            self.reference = image
            if np.isnan(image.any()):
                raise ValueError("Must have non-NaN image to " +\
                            "crosscorrelate for scale. Reference " +\
                            "image should also be centered.")

            self.refine_scaling(self.reference, centering=centering)
            self.pixscale_measured = self.pixscale_factor * self.pixel

            log.debug('measured pixel scale factor: %s', self.pixscale_factor)
            log.debug('measured pixel scale (mas): %s',\
                      leastsqnrm.rad2mas(self.pixscale_measured))

            self.fov = image.shape[0]

            self.fittingmodel = self.make_model(self.fov, bandpass=self.bandpass,
                                                 over=self.over, rotate=True,
                                                 centering=centering,
                                                 pixscale=self.pixscale_measured)

            log.debug('MODEL GENERATED AT PIXEL %s mas',
                      leastsqnrm.rad2mas(self.pixscale_measured))
            log.debug(' and ROTATION %s', self.rot_measured)

        else:
            self.fittingmodel = modelin

        self.soln, self.residual, self.cond =\
            leastsqnrm.matrix_operations(image, self.fittingmodel)

        log.debug('NRM_Model Raw Soln: %s', self.soln)

        self.rawDC = self.soln[-1]
        self.flux = self.soln[0]
        self.soln = self.soln / self.soln[0]
        self.fringeamp, self.fringephase =\
                        leastsqnrm.tan2visibilities(self.soln)
        self.piston = utils.fringes2pistons(self.fringephase, len(self.ctrs))
        self.closurephase = leastsqnrm.closurephase(self.fringephase, N=self.N)
        self.redundant_cps = leastsqnrm.redundant_cps(self.fringephase, N=self.N)
        self.redundant_cas = leastsqnrm.return_CAs(self.fringeamp, N=self.N)


    def create_modelpsf(self):
        """
        Short Summary
        -------------
        Make an image from the object's model and fit solutions, by setting the
        NRM_Model object's modelpsf attribute

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        try:
            self.modelpsf = np.zeros((self.fov, self.fov))
        except:
            self.modelpsf = np.zeros((self.fov_sim, self.fov_sim))

        for ind, coeff in enumerate(self.soln):
            self.modelpsf += self.flux * coeff * self.fittingmodel[:, :, ind]

        return None


    def refine_scaling(self, img, centering='PIXELCENTERED'):
        """
        Short Summary
        -------------
        Calculate pixel scale factor, rotation in addition to guess, and a
        goodness of fit parameter than can be compared in multiple iterations.
        Correlations are calculated in the image plane, in anticipation of
        data with many bad pixels.

        Parameters
        ----------
        img: 2D float array
            reference image

        centering: string
            type of centering; default='PIXELCENTERED'

        Returns
        -------
        Updated object's 'pixscale_factor': float
            Pixel scale factor

        Updated object's 'rot_measured': float
            Fit value of rotation

        Updated object's 'gof': float
            Goodness of fit value
        """

        if not hasattr(self, 'bandpass'):
            log.critical('There needs to be a specified bandpass/wavelength')

        reffov = img.shape[0]
        scallist = self.scallist
        scal_corrlist = np.zeros((len(scallist), reffov, reffov))
        pixscl_corrlist = scal_corrlist.copy()
        scal_corr = np.zeros(len(scallist))
        self.pixscl_corr = scal_corr.copy()

        self.pixscales = np.zeros(len(scallist))
        for q, scal in enumerate(scallist):

            self.test_pixscale = self.pixel * scal
            self.pixscales[q] = self.test_pixscale

            psf = self.simulate(bandpass=self.bandpass, fov=reffov,
                        pixel=self.test_pixscale, centering=centering)

            pixscl_corrlist[q, :, :] = run_data_correlate(img, psf)
            self.pixscl_corr[q] = np.max(pixscl_corrlist[q])
            if True in np.isnan(self.pixscl_corr):
                log.critical('Correlation produced NaNs,'\
                     ' please check your work!')

        self.pixscale_optimal, scal_maxy = utils.findmax(mag=self.pixscales,
                                vals=self.pixscl_corr)

        self.pixscale_factor = self.pixscale_optimal / self.pixel
        closestpixscale = self.pixscales[self.pixscl_corr ==\
                                          self.pixscl_corr.max()][0]

        log.debug('pixel scales search for max: %s %s ',
                  self.pixscales, self.pixscl_corr)
        log.debug('closestpixscale: %s ', closestpixscale)
        log.debug('NRM_Model pixel scale: %s', self.pixscale_optimal)
        log.debug('fraction of guess: %s', self.pixscale_factor)

        radlist = self.rotlist_rad
        corrlist = np.zeros((len(radlist), reffov, reffov))
        self.corrs = np.zeros(len(radlist))

        self.rots = self.rotlist_rad
        for q, rad in enumerate(radlist):
            psf = self.simulate(bandpass=self.bandpass, fov=reffov, \
                  pixel=self.pixscale_optimal, rotate=rad, \
                  centering=centering)

            corrlist[q, :, :] = run_data_correlate(psf, img)
            self.corrs[q] = np.max(corrlist[q])

        log.debug('rotation search for max: rots: %s corrs: %s', \
                            self.rots, self.corrs)

        self.rot_measured, maxy = utils.findmax(mag=self.rots, \
                                   vals=self.corrs)
        self.refpsf = self.simulate(bandpass=self.bandpass,
                                     pixel=self.pixscale_optimal,
                                     fov=reffov, rotate=self.rot_measured,
                                     centering=centering)

        log.info('------------ WHAT THE MATCHING ROUTINE FOUND ------------')
        log.info('scaling factor: %s', self.pixscale_factor)
        log.info('pixel scale (mas): %s', \
                        leastsqnrm.rad2mas(self.pixscale_optimal))
        log.info('rotation (rad): %s', self.rot_measured)
        log.info('--------------- -------------- ---------------')

        try:
            self.gof = goodness_of_fit(img, self.refpsf)
        except:
            log.warning("rewrite goodness_of_fit, it failed.")
            self.gof = False

        return self.pixscale_factor, self.rot_measured, self.gof


def makedisk(N, R, ctr=(0, 0), array="ODD"):
    """
    Short Summary
    -------------
    Calculate a 'disk', an array whose values =1 in a circular region near
    the center of the array, and =0 elsewhere. (Anand's emailed version)

    Parameters
    ----------
    N: integer

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

    if array == "ODD":
        M = (N - 1) / 2
        xx = np.linspace(-M - ctr[0], M - ctr[0], N)
        yy = np.linspace(-M - ctr[1], M - ctr[1], N)
    elif array == "EVEN":
        M = N / 2
        xx = np.linspace(-M - ctr[0], M - ctr[0] - 1, N)
        yy = np.linspace(-M - ctr[1], M - ctr[1] - 1, N)
    else: # unsupported value has been specified
        log.critical('An unsupported value has been specified in makedisk for array: %s', array)

    (x, y) = np.meshgrid(xx, yy.T)
    r = np.sqrt((x**2) + (y**2))
    array = np.zeros((N, N))
    array[r < R] = 1

    return array


def goodness_of_fit(data, bestfit, diskR=8):
    """
    Short Summary
    -------------
    Calculate goodness of fit between the data and the fit.

    Parameters
    ----------
    data: 2D float array
        input image ??

    bestfit: 2D float array
        fit to input image ??

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
        reference image ?

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
