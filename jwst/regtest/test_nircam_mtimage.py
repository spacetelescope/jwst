import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step


@pytest.mark.bigdata
@pytest.mark.parametrize("in_memory", [True, False])
def test_nircam_image_moving_target_i2d(rtdata, fitsdiff_default_kwargs, in_memory):
    """Test resampled i2d of moving target exposures for NIRCam imaging"""
    rtdata.get_asn("nircam/image/jw01252-o005_20240905t222322_image3_00001_asn.json")
    rtdata.output = "jw01252-o005_t003_nircam_clear-f277w_i2d.fits"
    args = ["calwebb_image3", rtdata.input, "--in_memory=" + str(in_memory)]
    Step.from_cmdline(args)
    rtdata.get_truth("truth/test_nircam_mtimage/jw01252-o005_t003_nircam_clear-f277w_i2d.fits")

    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 2e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
