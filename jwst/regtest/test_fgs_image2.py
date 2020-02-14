""" Test for the image2 pipeline for fgs data. This takes a rate file and
    generates flat field, i2d and cal files."""

import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.mark.bigdata
@pytest.mark.parametrize("output", [
    'jw86500007001_02101_00001_GUIDER2_cal.fits',
    'jw86500007001_02101_00001_GUIDER2_flat_field.fits',
    'jw86500007001_02101_00001_GUIDER2_i2d.fits'])
def test_fgs_image2(_jail, rtdata, fitsdiff_default_kwargs, output):
    """ Regression test for fgs data in the image2 pipeline"""
    rtdata.get_data("fgs/image2/jw86500007001_02101_00001_GUIDER2_rate.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_image2.cfg", rtdata.input,
            "--steps.flat_field.save_results=True",
            "--steps.resample.save_results=True"]
    Step.from_cmdline(args)
    rtdata.output = output

    rtdata.get_truth(os.path.join("truth/test_fgs_image2/", output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
