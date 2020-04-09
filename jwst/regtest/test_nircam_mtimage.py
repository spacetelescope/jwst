import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.mark.bigdata
def test_nircam_image_moving_target(rtdata, fitsdiff_default_kwargs):
    """Test resampled i2d of moving target exposures for NIRCam imaging"""
    collect_pipeline_cfgs("config")
    rtdata.get_asn("nircam/image/mt_asn.json")
    rtdata.output = "mt_assoc_i2d.fits"
    args = ["config/calwebb_image3.cfg", rtdata.input]
    Step.from_cmdline(args)
    rtdata.get_truth("truth/test_nircam_mtimage/mt_assoc_i2d.fits")

    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
