import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.mark.bigdata
def test_nirspec_image2_cal(_jail, rtdata, fitsdiff_default_kwargs):
    rtdata.get_data("nirspec/imaging/jw84600010001_02102_00001_nrs2_rate.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_image2.cfg", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "jw84600010001_02102_00001_nrs2_cal.fits"

    rtdata.get_truth("truth/test_nirspec_image2_cal/jw84600010001_02102_00001_nrs2_cal.fits")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
