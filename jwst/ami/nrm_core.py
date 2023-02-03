#! /usr/bin/env python

# Module for fringe fitting class
# Original python was by A. Greenbaum & A. Sivaramakrishnan

import os
import logging

import numpy as np
from scipy.special import comb
from scipy.stats import mstats

from stdatamodels.jwst import datamodels

from . import lg_model
from . import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class FringeFitter:
    def __init__(self, instrument_data, **kwargs):
        """
        Short Summary
        -------------
        For the given information on the instrument and mask, calculate the
        fringe observables (visibilities and closure phases in the image plane.

        Parameters
        ----------
        instrument_data - jwst.ami.instrument_data.NIRISS object
            Information on the mask geometry (namely # holes), instrument,
            wavelength obs mode.

        kwargs options:
            oversample - model oversampling (also how fine to measure the centering)
            psf_offset - subpixel centering of your data, if known
            npix - number of data pixels to use. Default is the shape of the data frame.
            find_rotation - will find the best pupil rotation that matches the data
        """
        self.instrument_data = instrument_data

        # Options
        if "oversample" in kwargs:
            self.oversample = kwargs["oversample"]
        else:
            # default oversampling is 3
            self.oversample = 3

        if "find_rotation" in kwargs:
            # can be True/False or 1/0
            self.find_rotation = kwargs["find_rotation"]
        else:
            self.find_rotation = False

        if "psf_offset_ff" in kwargs:  # if so do not find center of image in data
            self.psf_offset_ff = kwargs["psf_offset_ff"]
        else:
            self.psf_offset_ff = None  # find center of image in data

        if "npix" in kwargs:
            self.npix = kwargs["npix"]
        else:
            self.npix = 'default'

    def fit_fringes_all(self, input_model):
        """
        Short Summary
        ------------
        Extract the input data from input_model, and generate the best model to
        match the data (centering, scaling, rotation)

        Parameters
        ----------
        input_model: instance Data Model
            DM object for input

        Returns
        -------
        output_model: Fringe model object
            Fringe analysis data
        """

        self.scidata = self.instrument_data.read_data_model(input_model)

        nrm = lg_model.NrmModel(mask=self.instrument_data.mask,
                                pixscale=self.instrument_data.pscale_rad,
                                holeshape=self.instrument_data.holeshape,
                                affine2d=self.instrument_data.affine2d,
                                over=self.oversample)

        nrm.bandpass = self.instrument_data.wls[0]

        if self.npix == 'default':
            self.npix = self.scidata[:, :].shape[0]

        # New or modified in LG++
        # center the image on its peak pixel:
        # AS subtract 1 from "r" below  for testing >1/2 pixel offsets
        # AG 03-2019 -- is above comment still relevant?

        if self.instrument_data.arrname == "NIRC2_9NRM":
            self.ctrd = utils.center_imagepeak(self.scidata[0, :, :],
                                               r=(self.npix - 1) // 2 - 2, cntrimg=False)
        elif self.instrument_data.arrname == "gpi_g10s40":
            self.ctrd = utils.center_imagepeak(self.scidata[0, :, :],
                                               r=(self.npix - 1) // 2 - 2, cntrimg=True)
        else:
            self.ctrd = utils.center_imagepeak(self.scidata[:, :])

        nrm.reference = self.ctrd  # self. is the cropped image centered on the brightest pixel
        if self.psf_offset_ff is None:
            # returned values have offsets x-y flipped:
            # Finding centroids the Fourier way assumes no bad pixels case - Fourier domain mean slope
            centroid = utils.find_centroid(self.ctrd, self.instrument_data.threshold)  # offsets from brightest pixel ctr
            # use flipped centroids to update centroid of image for JWST - check parity for GPI, Vizier,...
            # pixel coordinates: - note the flip of [0] and [1] to match DS9 view

            nrm.xpos = centroid[1]  # flip 0 and 1 to convert
            nrm.ypos = centroid[0]  # flip 0 and 1
            nrm.psf_offset = nrm.xpos, nrm.ypos  # renamed .bestcenter to .psf_offset
        else:
            nrm.psf_offset = self.psf_offset_ff  # user-provided psf_offsetoffsets from array center are here.

        nrm.make_model(fov=self.ctrd.shape[0], bandpass=nrm.bandpass,
                       over=self.oversample,
                       psf_offset=nrm.psf_offset,
                       pixscale=nrm.pixel)

        nrm.fit_image(self.ctrd, modelin=nrm.model, psf_offset=nrm.psf_offset)

        """
        Attributes now stored in nrm object:

        -----------------------------------------------------------------------------
        soln            --- resulting sin/cos coefficients from least squares fitting
        fringephase     --- baseline phases in radians
        fringeamp       --- baseline amplitudes (flux normalized)
        redundant_cps   --- closure phases in radians
        redundant_cas   --- closure amplitudes
        residual        --- fit residuals [data - model solution]
        cond            --- matrix condition for inversion
        fringepistons   --- zero-mean piston opd in radians on each hole (eigenphases)
        -----------------------------------------------------------------------------
        """
        nrm.create_modelpsf()

        output_model = datamodels.AmiLgModel(
            fit_image=nrm.modelpsf,
            resid_image=nrm.residual,
            closure_amp_table=np.asarray(nrm.redundant_cas),
            closure_phase_table=np.asarray(nrm.redundant_cps),
            fringe_amp_table=np.asarray(nrm.fringeamp),
            fringe_phase_table=np.asarray(nrm.fringephase),
            pupil_phase_table=np.asarray(nrm.fringepistons),
            solns_table=np.asarray(nrm.soln))

        return output_model


class Calibrate:
    """
    A class to calibrate raw frame phases, closure phase, and squared visibilities.
    (The option to save to oifits files will be included in a future version.)
    """

    def __init__(self, objpaths, instrument_data, savedir=None, extra_dimension=None, **kwargs):
        """
        Summary
        -------
        Get statistics on each set of exposures, subtract calibrator phases from
        target phases, and divide target visibilities by cal visibilities.

        This will load all the observations into attributes:
        cp_mean_cal ... size [ncals, naxis2, ncp]
        cp_err_cal  ... size [ncals, naxis2, ncp]
        v2_mean_cal ... size [ncals, naxis2, nbl]
        v2_err_cal  ... size [ncals, naxis2, nbl]
        cp_mean_tar ... size [naxis2, nccp]
        cp_err_tar  ... size [naxis2, ncp]
        v2_mean_tar ... size [naxis2, nbl]
        v2_err_tar  ... size [naxis2, nbl]

        Parameters
        ----------
        objpaths: list (currently not used)
            list of directory paths (e.g. [tgtpth, calpth1, calpth2, ...]
            containing fringe observables for tgt, cal1, [cal2,...]. The
            first path is the target.  One or more calibrators follow.
            This parameter used to be the parameter 'paths'.

        instrument_data - jwst.ami.instrument_data.NIRISS object
            Information on the mask geometry (namely # holes), instrument,
            wavelength obs mode.

        savedir: string
            directory for output results (currently not used)

        extra_dimension: string (not currently used)
            dataset has an extra dimension for wavelength or polarization, for
            example. The value in each object level folder containing additional
            layer of data. Used to be parameter 'sub_dir_tag'.

        kwargs option:
            vflag: float
                visibility flag, to specify a cutoff and only accept exposures
                with higher visibilities
        """
        # Added 10/14/2016 -- here's a visibilities flag, where you can specify a
        # cutoff and only accept exposures with higher visibilities
        # Or maybe we should have it just reject outliers?
        # Let's have vflag be rejection of the fraction provided
        # e.g. if vflag=0.25, then exposures with the 25% lowest avg visibilities will be flagged
        if "vflag" in kwargs.keys():
            self.vflag = kwargs['vflag']
        else:
            self.vflag = 0.0

        # if no savedir specified, default is current working directory
        if savedir is None:
            self.savedir = os.getcwd()
        else:
            self.savedir = savedir

        try:
            os.listdir(savedir)
        except FileNotFoundError:
            os.mkdir(savedir)
        self.savedir = savedir

        # number of calibrators being used:
        self.ncals = len(objpaths) - 1  # number of calibrators, if zero, set to 1
        if self.ncals == 0:  # No calibrators given
            self.ncals = 1  # to avoid empty arrays
        self.nobjs = len(objpaths)  # number of total objects

        self.N = len(instrument_data.mask.ctrs)
        self.nbl = int(self.N * (self.N - 1) / 2)
        self.ncp = int(comb(self.N, 3))
        self.instrument_data = instrument_data

        # Additional axis (e.g., wavelength axis)
        # Can be size one. Defined by instrument_data wavelength array
        self.naxis2 = instrument_data.nwav

        # Set up all the arrays
        # Cal arrays have ncal axis, "wavelength" axis, and ncp axis
        self.cp_mean_cal = np.zeros((self.ncals, self.naxis2, self.ncp))
        self.cp_err_cal = np.zeros((self.ncals, self.naxis2, self.ncp))
        self.v2_mean_cal = np.zeros((self.ncals, self.naxis2, self.nbl))
        self.v2_err_cal = np.zeros((self.ncals, self.naxis2, self.nbl))
        self.pha_mean_cal = np.zeros((self.ncals, self.naxis2, self.nbl))
        self.pha_err_cal = np.zeros((self.ncals, self.naxis2, self.nbl))

        # target arrays have "wavelength" axis and ncp axis
        self.cp_mean_tar = np.zeros((self.naxis2, self.ncp))
        self.cp_err_tar = np.zeros((self.naxis2, self.ncp))
        self.v2_mean_tar = np.zeros((self.naxis2, self.nbl))
        self.v2_err_tar = np.zeros((self.naxis2, self.nbl))
        self.pha_mean_tar = np.zeros((self.naxis2, self.nbl))
        self.pha_err_tar = np.zeros((self.naxis2, self.nbl))

    # is there a subdirectory (e.g. for the exposure -- need to make this default)
        if extra_dimension is not None:
            self.extra_dimension = extra_dimension
            for ii in range(self.nobjs):
                exps = [f for f in os.listdir(objpaths[ii])
                        if self.extra_dimension in f and
                        os.path.isdir(os.path.join(objpaths[ii], f))]
                nexps = len(exps)

                amp = np.zeros((self.naxis2, nexps, self.nbl))
                pha = np.zeros((self.naxis2, nexps, self.nbl))
                cps = np.zeros((self.naxis2, nexps, self.ncp))

                # Create the cov matrix arrays
                if ii == 0:
                    self.cov_mat_tar = np.zeros((self.naxis2, nexps, nexps))
                    self.sigmasquared_tar = np.zeros((self.naxis2, nexps))
                    self.cov_mat_cal = np.zeros((self.naxis2, nexps, nexps))
                    self.sigmasquared_cal = np.zeros((self.naxis2, nexps))
                else:
                    pass
                expflag = []
                for qq in range(nexps):
                    # nwav files
                    cpfiles = [f for f in os.listdir(objpaths[ii] + exps[qq]) if "CPs" in f]

                    ampfiles = [f for f in os.listdir(objpaths[ii] + exps[qq])
                                if "amplitudes" in f]
                    phafiles = [f for f in os.listdir(objpaths[ii] + exps[qq]) if "phase" in f]

                    amp[0, qq, :] = np.loadtxt(objpaths[ii] + exps[qq] + "/" + ampfiles[0])
                    cps[0, qq, :] = np.loadtxt(objpaths[ii] + exps[qq] + "/" + cpfiles[0])
                    pha[0, qq, :] = np.loadtxt(objpaths[ii] + exps[qq] + "/" + phafiles[0])

                    # 10/14/2016 -- flag the exposure if we get amplitudes > 1
                    # Also flag the exposure if vflag is set, to reject fraction indicated
                    if True in (amp[:, qq, :] > 1):
                        expflag.append(qq)

                if self.vflag > 0.0:
                    self.ncut = int(self.vflag * nexps)  # how many are we cutting out
                    sorted_exps = np.argsort(amp.mean(axis=(0, -1)))
                    cut_exps = sorted_exps[:self.ncut]  # flag the ncut lowest exposures
                    expflag = expflag + list(cut_exps)

                # Create the cov matrix arrays
                if ii == 0:
                    rearr = np.rollaxis(cps, -1, 0).reshape(self.naxis2 * self.ncp, nexps)
                    R_i = rearr - rearr.mean(axis=1)[:, None]
                    R_j = R_i.T
                    self.cov = np.dot(R_i, R_j) / (nexps - 1)
                else:
                    rearr = np.rollaxis(cps, -1, 0).reshape(self.naxis2 * self.ncp, nexps)
                    R_i = rearr - rearr.mean(axis=1)[:, None]
                    R_j = R_i.T
                    self.cov += np.dot(R_i, R_j) / (nexps - 1)

                # Also adding a mask to calib steps

                if ii == 0:
                    # closure phases and squared visibilities
                    self.cp_mean_tar[0, :], self.cp_err_tar[0, :], \
                        self.v2_mean_tar[0, :], self.v2_err_tar[0, :], \
                        self.pha_mean_tar[0, :], self.pha_err_tar = \
                        self.calib_steps(cps[0, :, :], amp[0, :, :], pha[0, :, :], nexps, expflag=expflag)

                else:
                    # closure phases and visibilities
                    self.cp_mean_cal[ii - 1, 0, :], self.cp_err_cal[ii - 1, 0, :], \
                        self.v2_mean_cal[ii - 1, 0, :], self.v2_err_cal[ii - 1, 0, :], \
                        self.pha_mean_cal[ii - 1, 0, :], self.pha_err_cal[ii - 1, 0, :] = \
                        self.calib_steps(cps[0, :, :], amp[0, :, :], pha[0, :, :], nexps, expflag=expflag)
        else:
            for ii in range(self.nobjs):

                cpfiles = [f for f in os.listdir(objpaths[ii]) if "CPs" in f]
                ampfiles = [f for f in os.listdir(objpaths[ii]) if "amplitudes" in f]
                phafiles = [f for f in os.listdir(objpaths[ii]) if "phase" in f]
                nexps = len(cpfiles)

                amp = np.zeros((nexps, self.nbl))
                pha = np.zeros((nexps, self.nbl))
                cps = np.zeros((nexps, self.ncp))

                expflag = []
                for qq in range(nexps):
                    amp[qq, :] = np.loadtxt(objpaths[ii] + "/" + ampfiles[qq])
                    if True in (amp[qq, :] > 1):
                        expflag.append(qq)

                    pha[qq, :] = np.loadtxt(objpaths[ii] + "/" + phafiles[qq])
                    cps[qq, :] = np.loadtxt(objpaths[ii] + "/" + cpfiles[qq])

                # Covariance 06/27/2017
                if ii == 0:
                    rearr = np.rollaxis(cps, -1, 0).reshape(self.ncp, nexps)
                    R_i = rearr - rearr.mean(axis=1)[:, None]
                    R_j = R_i.T
                    self.cov = np.dot(R_i, R_j) / (nexps - 1)
                else:
                    rearr = np.rollaxis(cps, -1, 0).reshape(self.ncp, nexps)
                    R_i = rearr - rearr.mean(axis=1)[:, None]
                    R_j = R_i.T
                    self.cov += np.dot(R_i, R_j) / (nexps - 1)

                # Oct 14 2016 -- adding in a visibilities flag. Can't be >1 that doesn't make sense.
                # Also adding a mask to calib steps
                if ii == 0:
                    # closure phases and squared visibilities
                    self.cp_mean_tar[0, :], self.cp_err_tar[0, :], \
                        self.v2_mean_tar[0, :], self.v2_err_tar[0, :], \
                        self.pha_mean_tar[0, :], self.pha_err_tar[0, :] = \
                        self.calib_steps(cps, amp, pha, nexps, expflag=expflag)
                else:
                    # Fixed clunkiness!
                    # closure phases and visibilities
                    self.cp_mean_cal[ii - 1, 0, :], self.cp_err_cal[ii - 1, 0, :], \
                        self.v2_mean_cal[ii - 1, 0, :], self.v2_err_cal[ii - 1, 0, :], \
                        self.pha_mean_cal[ii - 1, 0, :], self.pha_err_cal[ii - 1, 0, :] = \
                        self.calib_steps(cps, amp, pha, nexps, expflag=expflag)

        # Combine mean calibrator values and errors
        self.cp_mean_tot = np.zeros(self.cp_mean_cal[0].shape)
        self.cp_err_tot = self.cp_mean_tot.copy()
        self.v2_mean_tot = np.zeros(self.v2_mean_cal[0].shape)
        self.v2_err_tot = self.v2_mean_tot.copy()
        self.pha_mean_tot = np.zeros(self.pha_mean_cal[0].shape)
        self.pha_err_tot = self.pha_mean_tot.copy()
        for ww in range(self.ncals):
            self.cp_mean_tot += self.cp_mean_cal[ww]
            self.cp_err_tot += self.cp_err_cal[ww]**2
            self.v2_mean_tot += self.v2_mean_cal[ww]
            self.v2_err_tot += self.v2_err_cal[ww]**2
            self.pha_mean_tot += self.pha_mean_cal[ww]
            self.pha_err_tot += self.pha_err_cal[ww]**2
        self.cp_mean_tot = self.cp_mean_tot / self.ncals
        self.cp_err_tot = np.sqrt(self.cp_err_tot)
        self.v2_mean_tot = self.v2_mean_tot / self.ncals
        self.v2_err_tot = np.sqrt(self.v2_err_tot)
        self.pha_mean_tot = self.pha_mean_tot / self.ncals
        self.pha_err_tot = np.sqrt(self.pha_err_tot)

        # Calibrate
        self.cp_calibrated = self.cp_mean_tar - self.cp_mean_tot
        self.cp_err_calibrated = np.sqrt(self.cp_err_tar**2 + self.cp_err_tot**2)
        self.v2_calibrated = self.v2_mean_tar / self.v2_mean_tot
        self.v2_err_calibrated = np.sqrt(self.v2_err_tar**2 + self.v2_err_tot**2)
        self.pha_calibrated = self.pha_mean_tar - self.pha_mean_tot
        self.pha_err_calibrated = np.sqrt(self.pha_err_tar**2 + self.pha_err_tot**2)

        # convert to degrees
        self.cp_calibrated_deg = self.cp_calibrated * 180 / np.pi
        self.cp_err_calibrated_deg = self.cp_err_calibrated * 180 / np.pi
        self.pha_calibrated_deg = self.pha_calibrated * 180 / np.pi
        self.pha_err_calibrated_deg = self.pha_err_calibrated * 180 / np.pi

    def calib_steps(self, cps, amps, pha, nexp, expflag=None):
        """
        Short Summary
        -------------
        Calculates closure phase and mean squared visibilities & standard error,
        and apply the exposure flag.

        Parameters
        ----------
        cps: 1D float array
            closure phases

        amps: 1D float array
            fringe visibility between each pair of holes

        pha: 1D float array
            phases

        nexp: integer
            number of exposures

        expflag: integer default=None
            number of flagged exposures

        Returns
        -------
        meancp: float
            mean of closure phases

        errcp: float
            error of closure phases

        meanv2: float
            mean squared visibilities

        errv2: float
            mean squared error of visibilities

        meanpha: float
            mean of phases

        errpha: float
            error of phases
        """

        # 10/14/16 Change flags exposures where vis > 1 anywhere
        # Apply the exposure flag
        expflag = None
        cpmask = np.zeros(cps.shape, dtype=bool)
        blmask = np.zeros(amps.shape, dtype=bool)
        if expflag is not None:
            nexp -= len(expflag)  # don't count the bad exposures
            cpmask[expflag, :] = True
            blmask[expflag, :] = True
        else:
            pass

        meancp = np.ma.masked_array(cps, mask=cpmask).mean(axis=0)
        meanv2 = np.ma.masked_array(amps, mask=blmask).mean(axis=0)**2
        meanpha = np.ma.masked_array(pha, mask=blmask).mean(axis=0)**2

        errcp = np.sqrt(mstats.moment(np.ma.masked_array(cps, mask=cpmask), moment=2, axis=0)) / np.sqrt(nexp)
        errv2 = np.sqrt(mstats.moment(np.ma.masked_array(amps**2, mask=blmask), moment=2, axis=0)) / np.sqrt(nexp)
        errpha = np.sqrt(mstats.moment(np.ma.masked_array(pha, mask=blmask), moment=2, axis=0)) / np.sqrt(nexp)

        # Set cutoff accd to Kraus 2008 - 2/3 of median
        errcp[errcp < (2 / 3.0) * np.median(errcp)] = (2 / 3.0) * np.median(errcp)
        errpha[errpha < (2 / 3.0) * np.median(errpha)] = (2 / 3.0) * np.median(errpha)
        errv2[errv2 < (2 / 3.0) * np.median(errv2)] = (2 / 3.0) * np.median(errv2)

        return meancp, errcp, meanv2, errv2, meanpha, errpha

    def save_to_txt(self):
        """
        Short Summary
        ------------
        Set up arrays to write to text files (not currently used)

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        tag = "_deg.txt"

        fns = ["cps" + tag, "cperr" + tag, "v2" + tag, "v2err" + tag]
        arrs = [self.cp_calibrated_deg, self.cp_err_calibrated_deg,
                self.v2_calibrated, self.v2_err_calibrated,
                self.pha_calibrated_deg, self.pha_err_calibrated_deg]
        self._save_txt(fns, arrs)

    def _save_txt(self, fns, arrays):
        """
        Short Summary
        -------------
        Write calibration arrays to text files (not currently used)

        Parameters
        ----------
        fns: list of strings
            names of 4 output files

        arrays: 1D float arrays
            4 arrays for closure phases, fringe visibilities, and errors

        Returns
        -------
        None
        """

        np.savetxt(self.savedir + "/" + fns[0], arrays[0])
        np.savetxt(self.savedir + "/" + fns[1], arrays[1])
        np.savetxt(self.savedir + "/" + fns[2], arrays[2])
        np.savetxt(self.savedir + "/" + fns[3], arrays[3])

        return None

    def save_to_oifits(self, **kwargs):
        """
        Short Summary
        -------------
        Save calibrated quantities to oifits files, with updated metadata
        (not yet supported).

        Parameters
        ----------
        kwargs options:
            phaseceil: float:
                value for keyword

            clip: integer ?
                clipping factor ?

            parang_range: float
                range of values of phase amplitudes ?

            avparang float
                average of range of values of phase amplitudes ?

        Returns
        -------
        None
        """

        # look for kwargs, e.g., phaseceil, anything else?
        if "phaseceil" in list(kwargs.keys()):
            self.phaseceil = kwargs["phaseceil"]
        else:
            # default for flagging closure phases (deg)
            self.phaseceil = 1.0e2  # degrees

        if "clip" in kwargs.keys():
            self.clip_wls = kwargs["clip"]
            nwav = self.naxis2 - 2 * self.clip_wls
            covshape = (self.naxis2 * self.ncp) - (2 * self.clip_wls * self.ncp)
            clippedcov = np.zeros((covshape, covshape))
            for k in range(self.ncp):
                clippedcov[nwav * k:nwav * (k + 1), nwav * k:nwav * (k + 1)] = \
                    self.cov[self.naxis2 * k + self.clip_wls:self.naxis2 * (k + 1) - self.clip_wls,
                             self.naxis2 * k + self.clip_wls:self.naxis2 * (k + 1) - self.clip_wls]
        else:
            # default is no clipping - maybe could set instrument-dependent clip in future
            self.clip_wls = None
            clippedcov = self.cov

        if not hasattr(self.instrument_data, "parang_range"):
            self.instrument_data.parang_range = 0.0
        if not hasattr(self.instrument_data, "avparang"):
            self.instrument_data.avparang = 0.0
        self.obskeywords = {
            'path': self.savedir + "/",
            'year': self.instrument_data.year,
            'month': self.instrument_data.month,
            'day': self.instrument_data.day,
            'TEL': self.instrument_data.telname,
            'instrument': self.instrument_data.instrument,
            'arrname': self.instrument_data.arrname,
            'object': self.instrument_data.objname,
            'RA': self.instrument_data.ra,
            'DEC': self.instrument_data.dec,
            'PARANG': self.instrument_data.avparang,
            'PARANGRANGE': self.instrument_data.parang_range,
            'PA': self.instrument_data.pa,
            'phaseceil': self.phaseceil,
            'covariance': clippedcov}
