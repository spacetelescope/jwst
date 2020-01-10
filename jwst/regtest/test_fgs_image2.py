""" Test for the image2 pipeline for fgs data. This takes a rate file and
    generates cal file."""

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.mark.bigdata
def test_fgs_image2(_jail, rtdata, fitsdiff_default_kwargs):
    rtdata.get_data("fgs/image/jw86500007001_02101_00001_GUIDER2_rate.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_image2.cfg", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "jw86500007001_02101_00001_GUIDER2_cal.fits"

    rtdata.get_truth("truth/test_fgs_image2/jw86500007001_02101_00001_GUIDER2_cal.fits")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
