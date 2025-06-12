"""Unit tests for AMI hextransformee module."""

import pytest
import numpy as np
from jwst.ami import instrument_data
from jwst.ami.utils import Affine2d


@pytest.mark.parametrize("chooseholes", (None, ["B2", "B4", "B5", "B6"]))
@pytest.mark.parametrize("affine2d", (None, Affine2d(rotradccw=0.4)))
@pytest.mark.parametrize("usebp", (True, False))
@pytest.mark.parametrize("firstfew", (None, 3))
@pytest.mark.parametrize("run_bpfix", (True, False))
def test_niriss(
    example_model,
    nrm_model_circular,
    bandpass,
    chooseholes,
    affine2d,
    usebp,
    firstfew,
    run_bpfix,
):
    """Test initialize NIRISS class with all available options."""

    filt = example_model.meta.instrument.filter
    niriss = instrument_data.NIRISS(
        filt,
        nrm_model_circular,
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

    assert scidata_ctrd.shape == dqmask_ctrd.shape
    if firstfew:
        assert scidata_ctrd.shape[0] == firstfew
        assert dqmask_ctrd.shape[0] == firstfew

    # test centering. should be centered around "real" source originally at (35, 35)
    shp_out = scidata_ctrd.shape[1:]
    assert shp_out == (63, 63)

    # Ensure real source is in the center of the image
    # Don't use zeroth integration because there's a bad pixel in that one, which
    # is not always masked out depending on input options
    assert np.unravel_index(np.argmax(scidata_ctrd[1]), shp_out) == (31, 31)

    # test bad pixels are fixed if run_bpfix is True
    if run_bpfix:
        assert np.max(scidata_ctrd) < 50  # bad pixel is fixed
    else:
        assert np.max(scidata_ctrd) == 100  # bad pixel retains original value of 100

    # test bad pixels are added to dq mask if usebp is True
    if run_bpfix and usebp:
        assert np.sum(dqmask_ctrd) == 1  # single bad pixel
    else:
        assert np.sum(dqmask_ctrd) == 0  # no masking

    # check attributes are set by read_data_model
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
