#
#  Module for defining data format, wavelength info, an mask geometry for these
#   instrument: NIRISS AMI
#

import logging
import numpy as np

from .mask_definitions import NRM_mask_definitions
from . import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

um = 1.0e-6


class NIRISS:
    def __init__(self, filt,
                 objname="obj",
                 src="A0V",
                 chooseholes=None,
                 affine2d=None,
                 bandpass=None,
                 **kwargs):
        """
        Short Summary
        ------------
        Initialize NIRISS class
        Either user has webbpsf and filter file can be read, or this will use a
          tophat and give a warning

        Parameters
        ----------
        filt: string
            filter name

        objname: string
            name of object

        src: string
            spectral type

        chooseholes: list
            None, or e.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask

        affine2d: Affine2d object
            Affine2d object

        bandpass: list
            None or [(wt,wlen),(wt,wlen),...].  Monochromatic would be e.g. [(1.0, 4.3e-6)]
            Explicit bandpass arg will replace *all* niriss filter-specific variables with
            the given bandpass, so you can simulate 21cm psfs through something called "F430M"!
        """
        self.chooseholes = chooseholes
        self.objname = objname
        self.filt = filt

        # 12 waves in f430 - data analysis:
        self.lam_bin = {"F277W": 50, "F380M": 20, "F430M": 40, "F480M": 30}

        # use 150 for 3 waves ax f430m; nominal values
        self.lam_c = {"F277W": 2.77e-6,  # central wavelength (SI)
                      "F380M": 3.8e-6,
                      "F430M": 4.28521033106325E-06,
                      "F480M": 4.8e-6}
        self.lam_w = {"F277W": 0.2, "F380M": 0.1, "F430M": 0.0436, "F480M": 0.08}  # fractional filter width

        self.throughput = utils.tophatfilter(self.lam_c[self.filt], self.lam_w[self.filt], npoints=11)

        # update nominal filter parameters with those of the filter read in and used in the analysis...
        # Weighted mean wavelength in meters, etc, etc "central wavelength" for the filter:
        from scipy.integrate import simps

        thru_st = np.stack(self.throughput, axis=1)
        thru_st_0 = thru_st[0, :]
        thru_st_1 = thru_st[1, :]

        num = (thru_st_0 * thru_st_1).sum()
        den = thru_st[0, :].sum()
        self.lam_c[self.filt] = num / den

        area = simps(thru_st_0, thru_st_1)
        ew = area / thru_st_0.max()  # equivalent width

        beta = ew / self.lam_c[self.filt]  # fractional bandpass
        self.lam_w[self.filt] = beta

        if bandpass is not None:
            bandpass = np.array(bandpass)  # type simplification
            wt = bandpass[:, 0]
            wl = bandpass[:, 1]
            cw = (wl * wt).sum() / wt.sum()  # Weighted mean wavelength in meters "central wavelength"
            area = simps(wt, wl)
            ew = area / wt.max()  # equivalent width
            beta = ew / cw  # fractional bandpass
            self.lam_c = {"F277W": cw, "F380M": cw, "F430M": cw, "F480M": cw, }
            self.lam_w = {"F277W": beta, "F380M": beta, "F430M": beta, "F480M": beta}
            self.throughput = bandpass

        self.wls = [self.throughput, ]
        # Wavelength info for NIRISS bands F277W, F380M, F430M, or F480M
        self.wavextension = ([self.lam_c[self.filt], ], [self.lam_w[self.filt], ])
        self.nwav = 1

        # only one NRM on JWST:
        self.telname = "NIRISS"
        self.instrument = "NIRISS"
        self.arrname = "jwst_g7s6c"
        self.holeshape = "hex"
        self.mask = NRM_mask_definitions(maskname=self.arrname, chooseholes=chooseholes, holeshape=self.holeshape)
        # save affine deformation of pupil object or create a no-deformation object.
        # We apply this when sampling the PSF, not to the pupil geometry.
        # This will set a default Ideal or a measured rotation, for example,
        # and include pixel scale changes due to pupil distortion.
        # Separating detector tilt pixel scale effects from pupil distortion effects is
        # yet to be determined... see comments in Affine class definition.
        if affine2d is None:
            self.affine2d = utils.Affine2d(mx=1.0, my=1.0,
                                           sx=0.0, sy=0.0,
                                           xo=0.0, yo=0.0, name="Ideal")
        else:
            self.affine2d = affine2d

        # finding centroid from phase slope only considered cv_phase data
        # when cv_abs data exceeds this cvsupport_threshold.
        # Absolute value of cv data normalized to unity maximum
        # for the threshold application.
        # Data reduction gurus: tweak the threshold value with experience...
        # Gurus: tweak cvsupport with use...
        self.cvsupport_threshold = {"F277W": 0.02, "F380M": 0.02, "F430M": 0.02, "F480M": 0.02}
        self.threshold = self.cvsupport_threshold[filt]

    def set_pscale(self, pscalex_deg=None, pscaley_deg=None):
        """
        Short Summary
        ------------
        Override pixel scale in header

        Parameters
        ----------
        pscalex_deg: float, degrees
            pixel scale in x-direction

        pscaley_deg: float, degrees
            pixel scale in y-direction

        Returns
        -------
        None

        """
        if pscalex_deg is not None:
            self.pscalex_deg = pscalex_deg
        if pscaley_deg is not None:
            self.pscaley_deg = pscaley_deg
        self.pscale_mas = 0.5 * (pscalex_deg + pscaley_deg) * (60 * 60 * 1000)
        self.pscale_rad = utils.mas2rad(self.pscale_mas)

    def read_data_model(self, input_model):
        """
        Short Summary
        -------------
        Retrieve info from input data model

        Parameters
        ----------
        input_model: instance Data Model
            DM object for input

        Returns
        -------
        Data and parameters from input data model
        """
        # The info4oif_dict will get pickled to disk when we write txt files of results.
        # That way we don't drag in objects like instrument_data into code that reads text results
        # and writes oifits files - a simple built-in dictionary is the only object used in this transfer.
        self.telname = "JWST"

        # To use ami_sim's eg 65.6 mas/pixel scale we hardcode it here.,,
        pscalex_deg = 65.6 / (1000 * 60 * 60)
        pscaley_deg = 65.6 / (1000 * 60 * 60)

        # Whatever we did set is averaged for isotropic pixel scale here
        self.pscale_mas = 0.5 * (pscalex_deg + pscaley_deg) * (60 * 60 * 1000)
        self.pscale_rad = utils.mas2rad(self.pscale_mas)
        self.mask = NRM_mask_definitions(maskname=self.arrname, chooseholes=self.chooseholes,
                                         holeshape=self.holeshape)

        return input_model.data

    def reset_nwav(self, nwav):
        """
        Reset self.nwav

        Parameters
        ----------
        nwav: integer
            length of axis3 for 3D input

        Returns
        -------
        """
        self.nwav = nwav
