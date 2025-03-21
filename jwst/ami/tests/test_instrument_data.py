"""Unit tests for AMI hextransformee module."""

import pytest

from jwst.ami import instrument_data
from jwst.ami.utils import Affine2d


@pytest.mark.parametrize("chooseholes", (None,)) #["B2", "B4", "B5", "B6"]))
@pytest.mark.parametrize("affine2d", (None, Affine2d(rotradccw=0.4)))
@pytest.mark.parametrize("usebp", (True, False))
@pytest.mark.parametrize("firstfew", (None, 3))
@pytest.mark.parametrize("run_bpfix", (True, False))
def test_niriss(example_model, nrm_model, bandpass, chooseholes, affine2d, usebp, firstfew, run_bpfix):
    """Test initialize NIRISS class with all available options."""

    filt = example_model.meta.instrument.filter 
    niriss = instrument_data.NIRISS(
        filt,
        nrm_model,
        chooseholes=chooseholes,
        affine2d=affine2d,
        bandpass=bandpass,
        usebp=usebp,
        firstfew=firstfew,
        run_bpfix=run_bpfix,
    )

    # check attributes are set
    # no nontrivial calculations in the init that are not already covered by other tests
    for att in [
        "filt",
        "nrm_model",
        "chooseholes",
        "affine2d",
        "throughput",
        "usebp",
        "firstfew",
        "run_bpfix",
        "lam_c",
        "lam_w",
        "nwav",
        "wls",
        "wavextension",
        "nwav",
        "telname",
        "instrument",
        "holeshape",
        "mask",
    ]:
        assert hasattr(niriss, att)

    scidata_ctrd, dqmask_ctrd = niriss.read_data_model(example_model)

    # TODOs:
    # test centering

    # test bad pixels are found and applied if run_bpfix is True

    # test bad pixel mask is returned as zeros if usebp is False,
    # and has actual DQ values if usebp is True

    # check attributes are set.
    # no more checks here for now, because this is likely to get refactored
    for att in [
        "pscale_rad",
        "pscale_mas",
        "roll_ref",
        "vparity",
        "v3iyang",
        "crpix1",
        "crpix2",
        "pupil",
        "proposer_name",
        "objname",
        "pi_name",
        "ra",
        "dec",
        "pmra",
        "pmdec",
        "ra_uncertainty",
        "dec_uncertainty",
        "date",
        "year",
        "month",
        "day",
        "itime",
        "ctrs_eqt",
        "ctrs_inst",
        "hdia",
        "nslices",
        "rootfn",
    ]:
        assert hasattr(niriss, att)