#! /usr/bin/env python
"""
RawOifits class: takes fringefitter class, which contains nrm_list and instrument_data attributes,
all info needed to write oifits.
populate structure needed to write out oifits files according to schema.
averaged and multi-integration versions, sigma-clipped stats over ints

CalibOifits class: takes two AmiOIModel datamodels and produces a final calibrated datamodel.
"""
import numpy as np

from scipy.special import comb
from astropy.stats import sigma_clipped_stats
from astropy.time.core import Time

from stdatamodels.jwst import datamodels
from . import leastsqnrm


class RawOifits:
    def __init__(self, fringefitter, method="mean"):
        """
        Class to store AMI data in the format required to write out to OIFITS files
        Angular quantities of input are in radians from fringe fitting; converted to degrees for saving.

        Parameters
        ----------
        fringefitter: FringeFitter object
            Object containing nrm_list attribute (list of nrm objects)
            and other info needed for OIFITS files
        method: string
            Method to average observables: mean or median. Default mean.

        Notes
        -----
        Based on ObservablesFromText from ImPlaneIA, e.g. 
        https://github.com/anand0xff/ImPlaneIA/blob/master/nrm_analysis/misctools/implane2oifits.py#L32
        """
        self.fringe_fitter = fringefitter
        self.n_holes = 7

        self.nslices = len(self.fringe_fitter.nrm_list)  # n ints
        self.n_baselines = int(comb(self.n_holes, 2))  # 21
        self.n_closure_phases = int(comb(self.n_holes, 3))  # 35
        self.n_closure_amplitudes = int(comb(self.n_holes, 4)) # also 35

        self.method = method

        self.ctrs_eqt = self.fringe_fitter.instrument_data.ctrs_eqt
        self.ctrs_inst = self.fringe_fitter.instrument_data.ctrs_inst
        self.pa = (
            self.fringe_fitter.instrument_data.pav3
        )  # header pav3, not including v3i_yang??

        self.bholes, self.bls = self._makebaselines()
        self.tholes, self.tuv = self._maketriples_all()
        self.qholes, self.quads = self._makequads_all()

    def make_obsarrays(self):
        """
        Make arrays of observables of the correct shape for saving to datamodels
        """
        # empty arrays of observables, (nslices,nobservables) shape.
        self.fringe_phases = np.zeros((self.nslices, self.n_baselines))
        self.fringe_amplitudes = np.zeros((self.nslices, self.n_baselines))
        self.closure_phases = np.zeros((self.nslices, self.n_closure_phases))
        self.t3_amplitudes = np.zeros((self.nslices, self.n_closure_phases))
        self.q4_phases = np.zeros((self.nslices, self.n_closure_amplitudes))
        self.closure_amplitudes = np.zeros((self.nslices, self.n_closure_amplitudes))
        self.pistons = np.zeros((self.nslices, self.n_holes))
        # model parameters
        self.solns = np.zeros((self.nslices, 44))

        # populate with each integration's observables
        for i, nrmslc in enumerate(self.fringe_fitter.nrm_list):
            self.fringe_phases[i, :] = nrmslc.fringephase  # FPs in radians
            self.fringe_amplitudes[i, :] = nrmslc.fringeamp
            self.closure_phases[i, :] = nrmslc.redundant_cps  # CPs in radians
            self.t3_amplitudes[i, :] = nrmslc.t3_amplitudes
            self.q4_phases[i, :] = nrmslc.q4_phases # quad phases in radians
            self.closure_amplitudes[i, :] = nrmslc.redundant_cas
            self.pistons[i, :] = nrmslc.fringepistons # segment pistons in radians
            self.solns[i, :] = nrmslc.soln

        self.fringe_amplitudes_squared = self.fringe_amplitudes ** 2  # squared visibilities

    def rotate_matrix(self, cov_mat, theta):
        """
        Rotate a covariance matrix by an angle.

        Parameters
        ----------
        cov_mat:
        theta: float
            Angle by which to rotate the matrix (radians)

        Returns
        -------

        """
        c, s = np.cos(theta), np.sin(theta)
        R_mat = [[c, -s],
                 [s, c]]
        # coordinate rotation from real/imaginary to absolute value/phase (modulus/argument)
        cv_rotated = np.linalg.multi_dot([np.transpose(R_mat), cov_mat, R_mat])
        return cv_rotated


    def average_observables(self, averfunc):
        """
        Calculate covariance matrices between fringe amplitudes/fringe phases, 
        and between triple product amps/closure phases, and closure amplitudes/quad phases.
        Convert r, theta (modulus, phase) to x,y. Calculate cov(x,y). Rotate resulting
        2x2 matrix back to r, theta. Take sqrt of relevant covariance element to be error.
        This must be done with phases in radians.

        Parameters
        ----------
        averfunc

        Returns
        -------
        avg_sqv, err_sqv, avg_fa, err_fa, avg_fp, err_fp, avg_cp, err_cp, avg_t3amp, err_t3amp, avg_ca, err_ca, avg_q4phi, err_q4phi, avg_pist, err_pist

        """
        covmats_fringes, covmats_triples, covmats_quads = self.observable_covariances(averfunc)

        if self.method == 'mean':
            avg_fa, _, std_fa = sigma_clipped_stats(self.fringe_amplitudes, axis=0)
            avg_fp, _, std_fp = sigma_clipped_stats(self.fringe_phases, axis=0)
            avg_sqv, _, std_sqv = sigma_clipped_stats(self.fringe_amplitudes**2, axis=0)
            avg_pist, _, err_pist = sigma_clipped_stats(self.pistons, axis=0)
        else:  # median
            _, avg_fa, std_fa = sigma_clipped_stats(self.fringe_amplitudes, axis=0)  # 21. std_fa is just for comparing to covariance
            _, avg_fp, std_fp  = sigma_clipped_stats(self.fringe_phases, axis=0)  # 21
            _, avg_sqv, std_sqv = sigma_clipped_stats(self.fringe_amplitudes**2, axis=0)
            _, avg_pist, err_pist = sigma_clipped_stats(self.pistons, axis=0)
        
        err_pist = err_pist/np.sqrt(self.nslices) # standard error of the mean
        err_fa, err_fp = self.err_from_covmat(covmats_fringes)

        # calculate squared visibility (fringe) amplitude uncertainties correctly
        err_sqv = (2 * avg_fa * err_fa) / np.sqrt(self.nslices)

        # calculate triple and quad quantities from **averaged** fringe amps and phases
        avg_t3amp = leastsqnrm.t3_amplitudes(avg_fa, n=self.n_holes)
        avg_cp = leastsqnrm.redundant_cps(avg_fp, n=self.n_holes)
        err_t3amp, err_cp = self.err_from_covmat(covmats_triples)
        
        avg_ca = leastsqnrm.closure_amplitudes(avg_fa, n=self.n_holes)
        avg_q4phi = leastsqnrm.q4_phases(avg_fp, n=self.n_holes)
        err_ca, err_q4phi = self.err_from_covmat(covmats_quads)

        return avg_sqv, err_sqv, avg_fa, err_fa, avg_fp, err_fp, avg_cp, err_cp, avg_t3amp, err_t3amp, avg_ca, err_ca, avg_q4phi, err_q4phi, avg_pist, err_pist

    def err_from_covmat(self, covmatlist):
        """
        Return sqrt of [0,0] and [1,1] elements of each of a list of covariance matrices,
        divided by sqrt(N_ints), for use as observable errors. (standard error of the mean)
        If using median, error calculation is questionable because this is NOT the standard
        error of the median
        
        Parameters
        ----------
        covmatlist: array
            array of covariance matrices for each baseline/triple/quad
            shape e.g. (21,2,2) or (35,2,2)

        Returns
        -------
        err_00: array
            standard errors of the mean of the first observable. shape e.g. (21)
        err_11: array
            standard errors of the mean of the second observable. shape e.g. (21)
        """
        err_00 = np.sqrt(np.array([covmat[0,0] for covmat in covmatlist]))/np.sqrt(self.nslices)
        err_11 = np.sqrt(np.array([covmat[1,1] for covmat in covmatlist]))/np.sqrt(self.nslices)

        return err_00, err_11

    def observable_covariances(self, averfunc):
        """
        input: nrm: ObservablesFromText instance
        Parameters
        ----------
        averfunc:

        Returns
        -------
        cov_mat_fringes:
        cov_mat_triples:
        cov_mat_quads:

        """
        # loop over 21 baselines
        cov_mat_fringes = []
        # these currently operate on all slices at once. other one operated on already-averaged slices
        for bl in np.arange(self.n_baselines):
            fringeamps = self.fringe_amplitudes[:,bl]
            fringephases = self.fringe_phases[:,bl]
            covmat = self.cov_r_theta(fringeamps, fringephases, averfunc)
            cov_mat_fringes.append(covmat)
        # loop over 35 triples
        cov_mat_triples = []
        for triple in np.arange(self.n_closure_phases):
            tripamp = self.t3_amplitudes[:,triple]
            triphase = self.closure_phases[:,triple]
            covmat = self.cov_r_theta(tripamp, triphase, averfunc)
            cov_mat_triples.append(covmat)
        # loop over 35 quads
        cov_mat_quads = []
        for quad in np.arange(self.n_closure_amplitudes):
            quadamp = self.closure_amplitudes[:,quad]
            quadphase = self.q4_phases[:,quad]
            covmat = self.cov_r_theta(quadamp, quadphase, averfunc)
            cov_mat_quads.append(covmat)

        # covmats to be written to oifits. store in rawoifits object? TBD
        # lists of cov mats have shape e.g. (21, 2, 2) or (35, 2, 2)

        return np.array(cov_mat_fringes), np.array(cov_mat_triples), np.array(cov_mat_quads)


    def cov_r_theta(self, rr, theta, averfunc):
        """
        Calculate covariance in x, y coordinates. Rotate covariance matrix by 
        **average** phase (over integrations) to get matrix in (r,theta)

        Parameters
        ----------
        rr: array
            complex number modulus 
        theta: array
            complex number phase
        averfunc: np.mean or np.median

        Returns
        -------
        cov_mat_r_theta: array (2,2)
            Covariance matrix in r, theta coordinates

        """
        xx = rr * np.cos(theta)
        yy = rr * np.sin(theta)
        cov_mat_xy = np.cov(xx, yy)
        cov_mat_r_theta = self.rotate_matrix(cov_mat_xy, averfunc(theta))
        return cov_mat_r_theta
        

    def make_oifits(self):
        """
        Perform final manipulations of observable arrays, calculate uncertainties, and
        populate AmiOIModel

        Returns
        -------
        m: AmiOIModel
            Fully populated datamodel

        """
        self.make_obsarrays()
        instrument_data = self.fringe_fitter.instrument_data
        observation_date = Time(
            "%s-%s-%s"
            % (instrument_data.year, instrument_data.month, instrument_data.day),
            format="fits",
        )

        # central wavelength, equiv. width from effstims used for fringe fitting
        wl = instrument_data.lam_c
        e_wl = instrument_data.lam_c * instrument_data.lam_w

        # Index 0 and 1 reversed to get the good u-v coverage (same fft)
        ucoord = self.bls[:, 1]
        vcoord = self.bls[:, 0]

        v1coord = self.tuv[:, 0, 0]
        u1coord = self.tuv[:, 0, 1]
        v2coord = self.tuv[:, 1, 0]
        u2coord = self.tuv[:, 1, 1]

        flagVis = [False] * self.n_baselines
        flagT3 = [False] * self.n_closure_phases

        # do the things done by populate_nrm here
        # average or don't, and get uncertainties
        # Unwrap phases
        shift2pi = np.zeros(self.closure_phases.shape)
        shift2pi[self.closure_phases >= 6] = 2 * np.pi
        shift2pi[self.closure_phases <= -6] = -2 * np.pi
        self.closure_phases -= shift2pi
        # this is problematic because phases were converted to DEGREES in make_obsarrays()?

        # Now we are setting up the observables to be written out to OIFITS
        if self.method not in ["mean", "median", "multi"]:
            self.method = "mean"  # default to mean
        # set these as attributes (some may exist and be overwritten)
        if self.method == "multi":
            self.vis2 = self.fringe_amplitudes_squared.T
            self.e_vis2 = np.zeros(self.vis2.shape)
            self.visamp = self.fringe_amplitudes.T
            self.e_visamp = np.zeros(self.visamp.shape)
            self.visphi = self.fringe_phases.T
            self.e_visphi = np.zeros(self.visphi.shape)
            self.closure_phases = self.closure_phases.T
            self.e_cp = np.zeros(self.closure_phases.shape)
            self.t3amp = self.t3_amplitudes.T
            self.e_t3amp = np.zeros(self.t3amp.shape)
            self.q4phi = self.q4_phases.T
            self.e_q4phi = np.zeros(self.q4_phases.shape)
            self.camp = self.closure_amplitudes.T
            self.e_camp = np.zeros(self.camp.shape)
            self.pist = self.pistons.T
            self.e_pist = np.zeros(self.pist.shape)


        elif self.method == "mean":
            self.vis2, self.e_vis2, self.visamp, self.e_visamp, self.visphi, self.e_visphi, self.closure_phases, self.e_cp, self.t3amp, self.e_t3amp, self.camp, self.e_camp, self.q4phi, self.e_q4phi, self.pist, self.e_pist =  self.average_observables(np.mean)

        else:  # take the median
            self.vis2, self.e_vis2, self.visamp, self.e_visamp, self.visphi, self.e_visphi, self.closure_phases, self.e_cp, self.t3amp, self.e_t3amp, self.camp, self.e_camp, self.q4phi, self.e_q4phi, self.pist, self.e_pist =  self.average_observables(np.median)
        
        # Convert angular quantities from radians to degrees
        self.visphi = np.rad2deg(self.visphi)
        self.e_visphi = np.rad2deg(self.e_visphi)
        self.closure_phases = np.rad2deg(self.closure_phases)
        self.e_cp = np.rad2deg(self.e_cp)
        self.q4phi = np.rad2deg(self.q4phi)
        self.e_q4phi = np.rad2deg(self.e_q4phi)
        self.pist = np.rad2deg(self.pist)
        self.e_pist = np.rad2deg(self.e_pist)


        # prepare arrays for OI_ARRAY ext
        self.staxy = instrument_data.ctrs_inst
        tel_name = ["A%i" % x for x in np.arange(self.n_holes) + 1]
        sta_name = tel_name
        diameter = [0] * self.n_holes

        staxyz = []
        for x in self.staxy:
            a = list(x)
            line = [a[0], a[1], 0]
            staxyz.append(line)

        sta_index = np.arange(self.n_holes) + 1

        pscale = instrument_data.pscale_mas / 1000.0  # arcsec
        isz = self.fringe_fitter.scidata.shape[
            1
        ]  # Size of the image to extract NRM data
        fov = [pscale * isz] * self.n_holes
        fovtype = ["RADIUS"] * self.n_holes

        m = datamodels.AmiOIModel()
        self.init_oimodel_arrays(m)

        # primary header keywords
        m.meta.telescope = instrument_data.telname
        m.meta.origin = "STScI"
        m.meta.instrument.name = instrument_data.instrument
        m.meta.program.pi_name = instrument_data.pi_name
        m.meta.target.proposer_name = instrument_data.proposer_name
        m.meta.observation.date = observation_date.fits
        m.meta.oifits.array_name = instrument_data.arrname
        m.meta.oifits.instrument_mode = instrument_data.pupil

        # oi_array extension data
        m.array["TEL_NAME"] = tel_name
        m.array["STA_NAME"] = sta_name
        m.array["STA_INDEX"] = sta_index
        m.array["DIAMETER"] = diameter
        m.array["STAXYZ"] = staxyz
        m.array["FOV"] = fov
        m.array["FOVTYPE"] = fovtype
        m.array["CTRS_EQT"] = instrument_data.ctrs_eqt
        m.array["PISTONS"] = self.pist
        m.array["PIST_ERR"] = self.e_pist

        # oi_target extension data
        m.target['TARGET_ID'] = [1]
        m.target['TARGET'] = instrument_data.objname
        m.target['RAEP0'] = instrument_data.ra
        m.target['DECEP0'] = instrument_data.dec
        m.target['EQUINOX'] = [2000]
        m.target['RA_ERR'] = instrument_data.ra_uncertainty
        m.target['DEC_ERR'] = instrument_data.dec_uncertainty
        m.target['SYSVEL'] = [0]
        m.target['VELTYP'] = ['UNKNOWN']
        m.target['VELDEF'] = ['OPTICAL']
        m.target['PMRA'] = instrument_data.pmra
        m.target['PMDEC'] = instrument_data.pmdec
        m.target['PMRA_ERR'] = [0]
        m.target['PMDEC_ERR'] = [0]
        m.target['PARALLAX'] = [0]
        m.target['PARA_ERR'] = [0]
        m.target['SPECTYP'] = ['UNKNOWN']

        # oi_vis extension data
        m.vis['TARGET_ID'] = 1
        m.vis['TIME'] = 0
        m.vis['MJD'] = observation_date.mjd
        m.vis['INT_TIME'] = instrument_data.itime
        m.vis['VISAMP'] = self.visamp
        m.vis['VISAMPERR'] = self.e_visamp
        m.vis['VISPHI'] = self.visphi
        m.vis['VISPHIERR'] = self.e_visphi
        m.vis['UCOORD'] = ucoord
        m.vis['VCOORD'] = vcoord
        m.vis['STA_INDEX'] = self._format_STAINDEX_V2(self.bholes)
        m.vis['FLAG'] = flagVis

        # oi_vis2 extension data
        m.vis2['TARGET_ID'] = 1
        m.vis2['TIME'] = 0
        m.vis2['MJD'] = observation_date.mjd
        m.vis2['INT_TIME'] = instrument_data.itime
        m.vis2['VIS2DATA'] = self.vis2
        m.vis2['VIS2ERR'] = self.e_vis2
        m.vis2['UCOORD'] = ucoord
        m.vis2['VCOORD'] = vcoord
        m.vis2['STA_INDEX'] = self._format_STAINDEX_V2(self.bholes)
        m.vis2['FLAG'] = flagVis

        # oi_t3 extension data
        m.t3['TARGET_ID'] = 1
        m.t3['TIME'] = 0
        m.t3['MJD'] = observation_date.mjd
        m.t3['T3AMP'] = self.t3amp
        m.t3['T3AMPERR'] = self.e_t3amp
        m.t3['T3PHI'] = self.closure_phases
        m.t3['T3PHIERR'] = self.e_cp
        m.t3['U1COORD'] = u1coord
        m.t3['V1COORD'] = v1coord
        m.t3['U2COORD'] = u2coord
        m.t3['V2COORD'] = v2coord
        m.t3['STA_INDEX'] = self._format_STAINDEX_T3(self.tholes)
        m.t3['FLAG'] = flagT3

        # oi_wavelength extension data
        m.wavelength["EFF_WAVE"] = wl
        m.wavelength["EFF_BAND"] = e_wl

        return m

    def init_oimodel_arrays(self, oimodel):
        """
        Set dtypes and initialize shapes for AmiOiModel arrays,
        depending on if averaged or multi-integration version.

        Parameters
        ----------
        oimodel: AmiOIModel object
            empty model
        """
        if self.method == "multi":
            # update dimensions of arrays for multi-integration oifits
            target_dtype = oimodel.target.dtype
            wavelength_dtype = np.dtype([("EFF_WAVE", "<f4"), ("EFF_BAND", "<f4")])
            array_dtype = np.dtype(
                [
                    ("TEL_NAME", "S16"),
                    ("STA_NAME", "S16"),
                    ("STA_INDEX", "<i2"),
                    ("DIAMETER", "<f4"),
                    ("STAXYZ", "<f8", (3,)),
                    ("FOV", "<f8"),
                    ("FOVTYPE", "S6"),
                    ("CTRS_EQT", "<f8", (2,)),
                    ("PISTONS", "<f8", (self.nslices,)),
                    ("PIST_ERR", "<f8", (self.nslices,)),
                ]
            )
            vis_dtype = np.dtype(
                [
                    ("TARGET_ID", "<i2"),
                    ("TIME", "<f8"),
                    ("MJD", "<f8"),
                    ("INT_TIME", "<f8"),
                    ("VISAMP", "<f8", (self.nslices,)),
                    ("VISAMPERR", "<f8", (self.nslices,)),
                    ("VISPHI", "<f8", (self.nslices,)),
                    ("VISPHIERR", "<f8", (self.nslices,)),
                    ("UCOORD", "<f8"),
                    ("VCOORD", "<f8"),
                    ("STA_INDEX", "<i2", (2,)),
                    ("FLAG", "i1"),
                ]
            )
            vis2_dtype = np.dtype(
                [
                    ("TARGET_ID", "<i2"),
                    ("TIME", "<f8"),
                    ("MJD", "<f8"),
                    ("INT_TIME", "<f8"),
                    ("VIS2DATA", "<f8", (self.nslices,)),
                    ("VIS2ERR", "<f8", (self.nslices,)),
                    ("UCOORD", "<f8"),
                    ("VCOORD", "<f8"),
                    ("STA_INDEX", "<i2", (2,)),
                    ("FLAG", "i1"),
                ]
            )
            t3_dtype = np.dtype(
                [
                    ("TARGET_ID", "<i2"),
                    ("TIME", "<f8"),
                    ("MJD", "<f8"),
                    ("INT_TIME", "<f8"),
                    ("T3AMP", "<f8", (self.nslices,)),
                    ("T3AMPERR", "<f8", (self.nslices,)),
                    ("T3PHI", "<f8", (self.nslices,)),
                    ("T3PHIERR", "<f8", (self.nslices,)),
                    ("U1COORD", "<f8"),
                    ("V1COORD", "<f8"),
                    ("U2COORD", "<f8"),
                    ("V2COORD", "<f8"),
                    ("STA_INDEX", "<i2", (3,)),
                    ("FLAG", "i1"),
                ]
            )
        else:
            target_dtype = oimodel.target.dtype
            wavelength_dtype = oimodel.wavelength.dtype
            array_dtype = np.dtype(
                [
                    ("TEL_NAME", "S16"),
                    ("STA_NAME", "S16"),
                    ("STA_INDEX", "<i2"),
                    ("DIAMETER", "<f4"),
                    ("STAXYZ", "<f8", (3,)),
                    ("FOV", "<f8"),
                    ("FOVTYPE", "S6"),
                    ("CTRS_EQT", "<f8", (2,)),
                    ("PISTONS", "<f8"),
                    ("PIST_ERR", "<f8"),
                ]
            )
            vis_dtype = oimodel.vis.dtype
            vis2_dtype = oimodel.vis2.dtype
            t3_dtype = oimodel.t3.dtype
        oimodel.array = np.zeros(self.n_holes, dtype=array_dtype)
        oimodel.target = np.zeros(1, dtype=target_dtype)
        oimodel.vis = np.zeros(self.n_baselines, dtype=vis_dtype)
        oimodel.vis2 = np.zeros(self.n_baselines, dtype=vis2_dtype)
        oimodel.t3 = np.zeros(self.n_closure_phases, dtype=t3_dtype)
        oimodel.wavelength = np.zeros(1, dtype=wavelength_dtype)

    def _maketriples_all(self):
        """
        Calculate all three-hole combinations, baselines

        Returns
        -------
        tarray: integer array
            Triple hole indices (0-indexed),
        float array of two uv vectors in all triangles
        """
        tlist = []
        uvlist = []
        for i in range(self.n_holes):
            for j in range(self.n_holes):
                for k in range(self.n_holes):
                    if i < j and j < k:
                        tlist.append((i, j, k))
                        uvlist.append(
                        (
                            self.ctrs_eqt[i] - self.ctrs_eqt[j],
                            self.ctrs_eqt[j] - self.ctrs_eqt[k],
                        )
                    )
        tarray = np.array(tlist).astype(int)
        return tarray, np.array(uvlist)

    def _makebaselines(self):
        """
        Calculate all hole pairs, baselines

        Returns
        -------
        barray: list
            Hole pairs indices, 0-indexed
        float array of baselines
        """
        blist = []
        bllist = []
        for i in range(self.n_holes):
            for j in range(self.n_holes):
                if i < j:
                    blist.append((i, j))
                    bllist.append(self.ctrs_eqt[i] - self.ctrs_eqt[j])
        return np.array(blist).astype(int), np.array(bllist)

    def _makequads_all(self):
        """ 
        Calculate all four-hole combinations (quads)
        
        Returns
        -------
        qarray: int array of four-hole quads (0-based)
        uvwlist: numpy array of u, v, w vectors for each quad
        """
        qlist = []
        uvwlist = []
        for i in range(self.n_holes):
            for j in range(self.n_holes):
                for k in range(self.n_holes):
                    for q in range(self.n_holes):
                        if i < j and j < k and k < q:
                            qlist.append((i, j, k, q))
                            uvwlist.append((self.ctrs_eqt[i] - self.ctrs_eqt[j],
                                            self.ctrs_eqt[j] - self.ctrs_eqt[k],
                                            self.ctrs_eqt[k] - self.ctrs_eqt[q]))
        qarray = np.array(qlist).astype(int)
        return qarray, np.array(uvwlist)

    def _format_STAINDEX_T3(self, tab):
        """
        Converts sta_index to save oifits T3 in the appropriate format

        Parameters
        ----------
        tab: array
            table of indices

        Returns
        -------
        sta_index: list
            Hole triples indices
        int array of triangles
        """    
        sta_index = []
        for x in tab:
            ap1 = int(x[0])
            ap2 = int(x[1])
            ap3 = int(x[2])
            if np.min(tab) == 0:
                line = np.array([ap1, ap2, ap3]) + 1
            else:
                line = np.array([ap1, ap2, ap3])
            sta_index.append(line)
        return sta_index

    def _format_STAINDEX_V2(self, tab):
        """
        Converts sta_index to save oifits V2 in the appropriate format

        Parameters
        ----------
        tab: array
            table of indices

        Returns
        -------
        sta_index: list
            Hole baseline indices
        int array of baselines
        """    
        sta_index = []
        for x in tab:
            ap1 = int(x[0])
            ap2 = int(x[1])
            if np.min(tab) == 0:
                line = np.array([ap1, ap2]) + 1 # RAC 2/2021
            else:
                line = np.array([ap1, ap2])
            sta_index.append(line)
        return sta_index


class CalibOifits:
    def __init__(self, targoimodel, caloimodel):
        """
        Calibrate (normalize) an AMI observation by subtracting closure phases
        of a reference star from those of a target and dividing visibility amplitudes
        of the target by those of the reference star.

        Parameters
        ----------
        targoimodel: AmiOIModlel, target
        caloimodel: AmiOIModlel, reference star (calibrator)
        """
        self.targoimodel = targoimodel
        self.caloimodel = caloimodel
        self.calib_oimodel = targoimodel.copy()

    def update_dtype(self):
        """
        Modify the dtype of OI array to include different pistons columns
        for calibrated OIFITS files
        """
        nrows = 7
        modified_dtype = np.dtype(
            [
                ("TEL_NAME", "S16"),
                ("STA_NAME", "S16"),
                ("STA_INDEX", "<i2"),
                ("DIAMETER", "<f4"),
                ("STAXYZ", "<f8", (3,)),
                ("FOV", "<f8"),
                ("FOVTYPE", "S6"),
                ("CTRS_EQT", "<f8", (2,)),
                ("PISTON_T", "<f8"),
                ("PISTON_C", "<f8"),
                ("PIST_ERR", "<f8"),
            ]
        )
        self.calib_oimodel.array = np.zeros(nrows, dtype=modified_dtype)

    def calibrate(self):
        """
        Apply the calibration (normalization) routine to calibrate the
        target AmiOIModel by the calibrator (reference star) AmiOIModel

        Returns
        -------
        calib_oimodel: AmiOIModel
            Calibrated AMI datamodel
        """
        cp_out = self.targoimodel.t3["T3PHI"] - self.caloimodel.t3["T3PHI"]
        sqv_out = self.targoimodel.vis2["VIS2DATA"] / self.caloimodel.vis2["VIS2DATA"]
        va_out = self.targoimodel.vis["VISAMP"] / self.caloimodel.vis["VISAMP"]
        # using standard propagation of error for multiplication/division
        # which assumes uncorrelated Gaussian errors (questionable)
        cperr_t = self.targoimodel.t3["T3PHIERR"]
        cperr_c = self.caloimodel.t3["T3PHIERR"]
        sqverr_c = self.targoimodel.vis2["VIS2ERR"]
        sqverr_t = self.caloimodel.vis2["VIS2ERR"]
        vaerr_t = self.targoimodel.vis["VISAMPERR"]
        vaerr_c = self.caloimodel.vis["VISAMPERR"]
        cperr_out = np.sqrt(cperr_t**2.0 + cperr_c**2.0)
        sqverr_out = sqv_out * np.sqrt(
            (sqverr_t / self.targoimodel.vis2["VIS2DATA"]) ** 2.0
            + (sqverr_c / self.caloimodel.vis2["VIS2DATA"]) ** 2.0
        )
        vaerr_out = va_out * np.sqrt(
            (vaerr_t / self.targoimodel.vis["VISAMP"]) ** 2.0
            + (vaerr_c / self.caloimodel.vis["VISAMP"]) ** 2.0
        )

        pistons_t = self.targoimodel.array["PISTONS"]
        pisterr_t = self.targoimodel.array["PIST_ERR"]
        pistons_c = self.caloimodel.array["PISTONS"]
        pisterr_c = self.caloimodel.array["PIST_ERR"]
        # sum in quadrature errors from target and calibrator pistons
        pisterr_out = np.sqrt(pisterr_t**2 + pisterr_c**2)

        # update OI array, which is currently all zeros, with input oi array
        # and updated piston columns
        self.update_dtype()
        self.calib_oimodel.array["TEL_NAME"] = self.targoimodel.array["TEL_NAME"]
        self.calib_oimodel.array["STA_NAME"] = self.targoimodel.array["STA_NAME"]
        self.calib_oimodel.array["STA_INDEX"] = self.targoimodel.array["STA_INDEX"]
        self.calib_oimodel.array["DIAMETER"] = self.targoimodel.array["DIAMETER"]
        self.calib_oimodel.array["STAXYZ"] = self.targoimodel.array["STAXYZ"]
        self.calib_oimodel.array["FOV"] = self.targoimodel.array["FOV"]
        self.calib_oimodel.array["FOVTYPE"] = self.targoimodel.array["FOVTYPE"]
        self.calib_oimodel.array["CTRS_EQT"] = self.targoimodel.array["CTRS_EQT"]
        self.calib_oimodel.array["PISTON_T"] = pistons_t
        self.calib_oimodel.array["PISTON_C"] = pistons_c
        self.calib_oimodel.array["PIST_ERR"] = pisterr_out

        # update calibrated oimodel arrays with calibrated observables
        self.calib_oimodel.t3["T3PHI"] = cp_out
        self.calib_oimodel.t3["T3PHIERR"] = cperr_out
        self.calib_oimodel.vis2["VIS2DATA"] = sqv_out
        self.calib_oimodel.vis2["VIS2ERR"] = sqverr_out
        self.calib_oimodel.vis["VISAMP"] = va_out
        self.calib_oimodel.vis["VISAMPERR"] = vaerr_out

        # add calibrated header keywords
        calname = self.caloimodel.meta.target.proposer_name  # name of calibrator star
        self.calib_oimodel.meta.oifits.calib = calname

        return self.calib_oimodel
