"""Unit tests for ami_normalize module and step."""

import pytest
import numpy as np
import stdatamodels.jwst.datamodels as dm

from jwst.ami import AmiNormalizeStep
from jwst.ami.utils import get_cw_beta
from .conftest import PXSC_MAS


RAW_AMP = 3
REF_AMP = 2
RAW_PHI = np.pi / 2
REF_PHI = -np.pi / 4
ERR = 0.1


@pytest.fixture
def oi_data(example_model, bandpass):
    """Create an AmiOiModel object with placeholder data."""
    oim = dm.AmiOIModel()
    # oim.update(example_model)

    nslices = 2
    n_holes = 7
    n_baselines = 21
    n_closure_phases = 35
    visamp = np.zeros((n_baselines)) + RAW_AMP
    visphi = np.zeros((n_baselines)) + RAW_PHI
    t3amp = np.zeros((n_closure_phases)) + RAW_AMP
    t3phi = np.zeros((n_closure_phases)) + RAW_PHI
    sta_index = np.arange(n_holes) + 1
    pistons = np.zeros((n_holes,))
    flag_vis = [False] * n_baselines
    flag_t3 = [False] * n_closure_phases

    isz = example_model.data.shape[1]
    pscale = PXSC_MAS / 1000.0  # arcsec
    lam_c, lam_w = get_cw_beta(bandpass)

    # oi_header data
    oim.meta.oifits.array_name = "jwst_ami"
    oim.meta.oifits.instrument_mode = example_model.meta.instrument.pupil

    # oimodel.array = np.zeros(self.n_holes, dtype=array_dtype)
    # oimodel.target = np.zeros(1, dtype=target_dtype)
    # oimodel.wavelength = np.zeros(1, dtype=wavelength_dtype)

    # oi_array extension data
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
    oim.array = np.zeros(n_holes, dtype=array_dtype)
    oim.array["TEL_NAME"] = [f"A{x:d}" for x in sta_index]
    oim.array["STA_NAME"] = [f"A{x:d}" for x in sta_index]
    oim.array["STA_INDEX"] = sta_index
    oim.array["DIAMETER"] = [0] * n_holes
    oim.array["STAXYZ"] = [
        [0, 0, 0],
    ] * n_holes
    oim.array["FOV"] = [pscale * isz] * n_holes
    oim.array["FOVTYPE"] = ["RADIUS"] * n_holes
    # oim.array["CTRS_EQT"] = instrument_data.ctrs_eqt
    oim.array["PISTONS"] = pistons
    oim.array["PIST_ERR"] = pistons * ERR

    # oi_target extension data
    oim.target["TARGET_ID"] = [1]
    oim.target["TARGET"] = "target"
    oim.target["RAEP0"] = example_model.meta.target.ra
    oim.target["DECEP0"] = example_model.meta.target.dec
    oim.target["EQUINOX"] = [2000]
    oim.target["RA_ERR"] = example_model.meta.target.ra_uncertainty
    oim.target["DEC_ERR"] = example_model.meta.target.dec_uncertainty
    oim.target["SYSVEL"] = [0]
    oim.target["VELTYP"] = ["UNKNOWN"]
    oim.target["VELDEF"] = ["OPTICAL"]
    oim.target["PMRA"] = example_model.meta.target.proper_motion_ra
    oim.target["PMDEC"] = example_model.meta.target.proper_motion_dec
    oim.target["PMRA_ERR"] = [0]
    oim.target["PMDEC_ERR"] = [0]
    oim.target["PARALLAX"] = [0]
    oim.target["PARA_ERR"] = [0]
    oim.target["SPECTYP"] = ["UNKNOWN"]

    # oi_vis extension data
    oim.vis = np.zeros(n_baselines, dtype=oim.vis.dtype)
    oim.vis["TARGET_ID"] = 1
    oim.vis["TIME"] = 0
    # oim.vis["MJD"] = 0
    # oim.vis["INT_TIME"] = 0
    oim.vis["VISAMP"] = visamp
    oim.vis["VISAMPERR"] = visamp * ERR
    oim.vis["VISPHI"] = visphi
    oim.vis["VISPHIERR"] = visphi * ERR
    # oim.vis["UCOORD"] = ucoord
    # oim.vis["VCOORD"] = vcoord
    # oim.vis["STA_INDEX"] = sta_index
    oim.vis["FLAG"] = flag_vis

    # oi_vis2 extension data
    oim.vis2 = np.zeros(n_baselines, dtype=oim.vis2.dtype)
    oim.vis2["TARGET_ID"] = 1
    oim.vis2["TIME"] = 0
    # oim.vis2["MJD"] = observation_date.mjd
    # oim.vis2["INT_TIME"] = instrument_data.itime
    oim.vis2["VIS2DATA"] = (visamp**2).T
    oim.vis2["VIS2ERR"] = (visamp**2).T * ERR
    # oim.vis2["UCOORD"] = ucoord
    # oim.vis2["VCOORD"] = vcoord
    # oim.vis2["STA_INDEX"] = np.arange(n_holes) + 1
    oim.vis2["FLAG"] = flag_vis

    # oi_t3 extension data
    oim.t3 = np.zeros(n_closure_phases, dtype=oim.t3.dtype)
    oim.t3["TARGET_ID"] = 1
    oim.t3["TIME"] = 0
    # oim.t3["MJD"] = observation_date.mjd
    oim.t3["T3AMP"] = t3amp
    oim.t3["T3AMPERR"] = t3amp * ERR
    oim.t3["T3PHI"] = t3phi
    oim.t3["T3PHIERR"] = t3phi * ERR
    # oim.t3["U1COORD"] = u1coord
    # oim.t3["V1COORD"] = v1coord
    # oim.t3["U2COORD"] = u2coord
    # oim.t3["V2COORD"] = v2coord
    # oim.t3["STA_INDEX"] = sta_index
    oim.t3["FLAG"] = flag_t3

    # oi_wavelength extension data
    oim.wavelength["EFF_WAVE"] = lam_c
    oim.wavelength["EFF_BAND"] = lam_c * lam_w

    return oim


@pytest.fixture
def ref_data(oi_data):
    ref_data = oi_data.copy()
    ref_data.vis["VISAMP"] = np.zeros_like(oi_data.vis["VISAMP"]) + REF_AMP
    ref_data.vis["VISPHI"] = np.zeros_like(oi_data.vis["VISPHI"]) + REF_PHI
    ref_data.vis2["VIS2DATA"] = np.zeros_like(oi_data.vis2["VIS2DATA"]) + REF_AMP**2
    ref_data.t3["T3AMP"] = np.zeros_like(oi_data.t3["T3AMP"]) + REF_AMP
    ref_data.t3["T3PHI"] = np.zeros_like(oi_data.t3["T3PHI"]) + REF_PHI

    return ref_data


def test_ami_normalize(oi_data, ref_data):
    """
    Test the AmiNormalizeStep.

    This test provides coverage for AmiNormalizeStep, ami_normalize.py,
    and the CalibOiFits class inside oifits.py, because internally it
    initializes a CalibOiFits class and then calls its calibrate() method.
    """

    result = AmiNormalizeStep.call(oi_data, ref_data)

    assert isinstance(result, dm.AmiOIModel)
    assert np.allclose(result.vis["VISAMP"], RAW_AMP / REF_AMP)
    assert np.allclose(result.vis["VISPHI"], RAW_PHI)
    assert np.allclose(result.vis2["VIS2DATA"], (RAW_AMP / REF_AMP) ** 2)
    assert np.allclose(result.t3["T3AMP"], RAW_AMP)
    assert np.allclose(result.t3["T3PHI"], RAW_PHI - REF_PHI)
