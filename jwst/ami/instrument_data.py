#
#  Module for defining data format, wavelength info, an mask geometry for these
#   instrument: NIRISS AMI
#

import logging
import numpy as np

from .mask_definitions import NRM_mask_definitions
from . import utils
from . import bp_fix
from stdatamodels.jwst.datamodels import dqflags
import copy


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class NIRISS:
    def __init__(self,
                 filt,
                 nrm_model,
                 chooseholes=None,
                 affine2d=None,
                 bandpass=None,
                 usebp=True,
                 firstfew=None,
                 rotsearch_parameters=None,
                 oversample=None,
                 psf_offset=None,
                 run_bpfix=True
                 ):
        """
        Initialize NIRISS class

        Parameters
        ----------
        filt: string
            filter name

        chooseholes: list
            None, or e.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask

        affine2d: Affine2d object
            Affine2d object

        bandpass: synphot spectrum or array
            None, synphot object or [(wt,wlen),(wt,wlen),...].  Monochromatic would be e.g. [(1.0, 4.3e-6)]
            Explicit bandpass arg will replace *all* niriss filter-specific variables with
            the given bandpass, so you could simulate, for example, a 21cm psf through something called "F430M"!

        usebp : boolean
            If True, exclude pixels marked DO_NOT_USE from fringe fitting

        firstfew : integer
            If not None, process only the first few integrations

        chooseholes : string
            If not None, fit only certain fringes e.g. ['B4','B5','B6','C2']

        affine2d : Affine2D object
            None or user-defined Affine2d object

        run_bpfix : boolean
            Run Fourier bad pixel fix on cropped data
        """
        self.run_bpfix = run_bpfix
        self.usebp = usebp
        self.chooseholes = chooseholes
        self.filt = filt
        self.throughput = bandpass
        self.firstfew = firstfew
        self.nrm_model = nrm_model

        self.lam_c, self.lam_w = utils.get_cw_beta(self.throughput)
        self.wls = [
            self.throughput,
        ]
        # Wavelength info for NIRISS bands F277W, F380M, F430M, or F480M
        self.wavextension = (
            [
                self.lam_c,
            ],
            [
                self.lam_w,
            ],
        )
        self.nwav = 1

        # only one NRM on JWST:
        self.telname = "JWST"
        self.instrument = "NIRISS"
        self.arrname = "jwst_g7s6c"
        self.holeshape = "hex"
        self.mask = NRM_mask_definitions(
            maskname=self.arrname,
            chooseholes=self.chooseholes,
            holeshape=self.holeshape,
        )

        # save affine deformation of pupil object or create a no-deformation object.
        # We apply this when sampling the PSF, not to the pupil geometry.
        # This will set a default Ideal or a measured rotation, for example,
        # and include pixel scale changes due to pupil distortion.
        # Separating detector tilt pixel scale effects from pupil distortion effects is
        # yet to be determined... see comments in Affine class definition.
        if affine2d is None:
            self.affine2d = utils.Affine2d(
                mx=1.0, my=1.0, sx=0.0, sy=0.0, xo=0.0, yo=0.0, name="Ideal"
            )
        else:
            self.affine2d = affine2d

        # finding centroid from phase slope only considered cv_phase data
        # when cv_abs data exceeds this cvsupport_threshold.
        # Absolute value of cv data normalized to unity maximum
        # for the threshold application.
        # Data reduction gurus: tweak the threshold value with experience...
        # Gurus: tweak cvsupport with use...
        self.cvsupport_threshold = {
            "F277W": 0.02,
            "F380M": 0.02,
            "F430M": 0.02,
            "F480M": 0.02,
        }
        self.threshold = self.cvsupport_threshold[filt]

    def set_pscale(self, pscalex_deg=None, pscaley_deg=None):
        """
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
        Retrieve info from input data model and store in NIRISS class.
        Trim refpix and roughly center science data and dq array.
        Run Fourier bad pixel correction before returning science data.

        Parameters
        ----------
        input_model: instance Data Model
            DM object for input

        Returns
        -------
        scidata_ctrd: numpy array
            Cropped, centered, optionally cleaned AMI data
        dqmask_ctrd:
            Cropped, centered mask of bad pixels
        """

        # all instrumentdata attributes will be available when oifits files written out
        scidata = copy.deepcopy(np.array(input_model.data))
        bpdata = copy.deepcopy(np.array(input_model.dq))

        # pixel scale recalculated and averaged
        pscaledegx, pscaledegy = utils.degrees_per_pixel(input_model)
        # At some point we want to use different X and Y pixel scales
        pscale_deg = np.mean([pscaledegx, pscaledegy])
        self.pscale_rad = np.deg2rad(pscale_deg)
        self.pscale_mas = pscale_deg * (60 * 60 * 1000)
        self.pav3 = input_model.meta.pointing.pa_v3
        self.vparity = input_model.meta.wcsinfo.vparity
        self.v3iyang = input_model.meta.wcsinfo.v3yangle
        self.parangh = input_model.meta.wcsinfo.roll_ref
        self.crpix1 = input_model.meta.wcsinfo.crpix1
        self.crpix2 = input_model.meta.wcsinfo.crpix2
        self.pupil = input_model.meta.instrument.pupil
        self.proposer_name = input_model.meta.target.proposer_name
        if input_model.meta.target.catalog_name == "UNKNOWN":
            objname = input_model.meta.target.proposer_name
        else:
            objname = input_model.meta.target.catalog_name
        self.objname = objname.replace("-", " ")
        self.pi_name = input_model.meta.program.pi_name
        self.ra = input_model.meta.target.ra
        self.dec = input_model.meta.target.dec
        self.pmra = input_model.meta.target.proper_motion_ra
        self.pmdec = input_model.meta.target.proper_motion_dec
        self.ra_uncertainty = input_model.meta.target.ra_uncertainty
        self.dec_uncertainty = input_model.meta.target.dec_uncertainty


        datestr = input_model.meta.visit.start_time.replace(" ", "T")
        self.date = datestr  # is this the right start time?
        self.year = datestr[:4]
        self.month = datestr[5:7]
        self.day = datestr[8:10]
        effinttm = input_model.meta.exposure.effective_exposure_time
        nints = input_model.meta.exposure.nints
        # if 2d input, model has already been expanded to 3d, so check 0th dimension
        if input_model.data.shape[0] == 1:
            self.itime = effinttm * nints
        else:
            self.itime = effinttm
            if self.firstfew is not None:
                if scidata.shape[0] > self.firstfew:
                    log.info(f"Analyzing only the first {self.firstfew:d} integrations")
                    scidata = scidata[: self.firstfew, :, :]
                    bpdata = bpdata[: self.firstfew, :, :]
                else:
                    log.warning(f"Input firstfew={self.firstfew:d} is greater than the number of integrations")
                    log.warning("All integrations will be analyzed")
            self.nwav = scidata.shape[0]
            [self.wls.append(self.wls[0]) for f in range(self.nwav - 1)]
        # Rotate mask hole centers by pav3 + v3i_yang to be in equatorial coordinates
        ctrs_sky = self.mast2sky()
        oifctrs = np.zeros(self.mask.ctrs.shape)
        oifctrs[:, 0] = ctrs_sky[:, 1].copy() * -1
        oifctrs[:, 1] = ctrs_sky[:, 0].copy() * -1
        self.ctrs_eqt = oifctrs
        self.ctrs_inst = self.mask.ctrs
        self.hdia = self.mask.hdia
        self.nslices = self.nwav

        # Trim refpix from all slices
        scidata = scidata[:, 4:, :]
        bpdata = bpdata[:, 4:, :]

        # find peak in median of refpix-trimmed scidata
        med_im = np.median(scidata, axis=0)

        # Use median image to find big CR hits not already caught by pipeline
        std_im = np.std(scidata, axis=0)
        mediandiff = np.empty_like(scidata)
        mediandiff[:, :, :] = scidata - med_im
        nsigma = 10
        outliers = np.where(mediandiff > nsigma * std_im)
        outliers2 = np.argwhere(mediandiff > nsigma * std_im)

        dqvalues = bpdata[outliers]
        log.info(
            f"{len(dqvalues)} additional pixels >10-sig from median of stack found"
        )
        # decompose DQ values to check if they are already flagged DNU
        count = 0
        for loc, dq_value in zip(outliers2, dqvalues):
            bitarr = np.binary_repr(dq_value)
            bad_types = []
            for i, elem in enumerate(bitarr[::-1]):
                if elem == str(1):
                    badval = 2**i
                    key = next(
                        key for key, value in dqflags.pixel.items() if value == badval
                    )
                    bad_types.append(key)
            if "DO_NOT_USE" not in bad_types:
                bpdata[loc[0], loc[1], loc[2]] += 1
                count += 1
        log.info(f"{count} DO_NOT_USE flags added to DQ array for found outlier pixels")

        # Roughly center scidata, bpdata around peak pixel position
        peakx, peaky, r = utils.min_distance_to_edge(med_im)
        scidata_ctrd = scidata[
            :, int(peakx - r):int(peakx + r + 1), int(peaky - r):int(peaky + r + 1)
        ]
        bpdata_ctrd = bpdata[
            :, int(peakx - r):int(peakx + r + 1), int(peaky - r):int(peaky + r + 1)
        ]

        log.info(
            "Cropping all integrations to %ix%i pixels around peak (%i,%i)"
            % (2 * r + 1, 2 * r + 1, peakx + 4, peaky)
        )  # +4 because of trimmed refpx
        # apply bp fix here
        if self.run_bpfix:
            log.info('Applying Fourier bad pixel correction to cropped data, updating DQ array')
            scidata_ctrd, bpdata_ctrd = bp_fix.fix_bad_pixels(
                scidata_ctrd,
                bpdata_ctrd,
                input_model.meta.instrument.filter,
                self.pscale_mas,
                self.nrm_model,
            )
        else:
            log.info("Not running Fourier bad pixel fix")

        self.rootfn = input_model.meta.filename.replace(".fits", "")

        # all info needed to write out oifits should be stored in NIRISS object attributes

        # Make a bad pixel mask, either from real DQ data or zeros if usebp=False
        if self.usebp:
            log.info(
                "usebp flag set to TRUE: bad pixels will be excluded from model fit"
            )
            DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
            JUMP_DET = dqflags.pixel["JUMP_DET"]
            dq_dnu = bpdata_ctrd & DO_NOT_USE == DO_NOT_USE
            dq_jump = bpdata_ctrd & JUMP_DET == JUMP_DET
            dqmask_ctrd = dq_dnu | dq_jump
        else:
            log.info("usebp flag set to FALSE: all pixels will be used in model fit")
            dqmask_ctrd = np.zeros_like(scidata_ctrd)

        return scidata_ctrd, dqmask_ctrd

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

    def mast2sky(self):
        """
        Rotate hole center coordinates:
            Clockwise by the V3 position angle - V3I_YANG from north in degrees if VPARITY = -1
            Counterclockwise by the V3 position angle - V3I_YANG from north in degrees if VPARITY = 1
        Hole center coords are in the V2, V3 plane in meters.

        Returns
        -------
        ctrs_rot: array
            Rotated coordinates to be put in OIFITS files.
        """
        mask_ctrs = copy.deepcopy(self.mask.ctrs)
        # rotate by an extra 90 degrees (RAC 9/21)
        # these coords are just used to orient output in OIFITS files
        # NOT used for the fringe fitting itself
        mask_ctrs = utils.rotate2dccw(mask_ctrs, np.pi / 2.0)
        vpar = self.vparity  # Relative sense of rotation between Ideal xy and V2V3
        rot_ang = self.pav3 - self.v3iyang  # subject to change!

        if self.pav3 == 0.0:
            return mask_ctrs
        else:
            # Using rotate2sccw, which rotates **vectors** CCW in a fixed coordinate system,
            # so to rotate coord system CW instead of the vector, reverse sign of rotation angle.  Double-check comment
            if vpar == -1:
                # rotate clockwise  <rotate coords clockwise>
                ctrs_rot = utils.rotate2dccw(mask_ctrs, np.deg2rad(-rot_ang))
                log.info(
                    f"Rotating mask hole centers clockwise by {rot_ang:.3f} degrees"
                )
            else:
                # counterclockwise  <rotate coords counterclockwise>
                ctrs_rot = utils.rotate2dccw(mask_ctrs, np.deg2rad(rot_ang))
                log.info(
                    f"Rotating mask hole centers counterclockwise by {rot_ang:.3f} degrees"
                )
            return ctrs_rot
