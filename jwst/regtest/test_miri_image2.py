import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.mark.bigdata
def test_miri_image2_cal(_jail, rtdata, fitsdiff_default_kwargs):
    rtdata.get_data("miri/image/jw00001001001_01101_00001_mirimage_rate.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_image2.cfg", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "jw00001001001_01101_00001_mirimage_cal.fits"

    rtdata.get_truth("truth/test_miri_image2_cal/jw00001001001_01101_00001_mirimage_cal.fits")

    fitsdiff_default_kwargs["rtol"] = 0.0001
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_image2 pipeline on MIRI imaging data."""
    rtdata = rtdata_module
    rtdata.get_data("miri/image/jw00001001001_01101_00001_mirimage_rate.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_image2.cfg", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("output",[
    "jw00001001001_01101_00001_mirimage_cal.fits",
    "jw00001001001_01101_00001_mirimage_i2d.fits"],
    ids=["cal", "i2d"])
def test_miri_image2(run_pipeline, fitsdiff_default_kwargs, output):
    """Regression test of calwebb_image2 pipeline performed on MIRI data."""
    rtdata = run_pipeline
    rtdata.output = output

    rtdata.get_truth(os.path.join("truth/test_miri_image2_cal", output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
