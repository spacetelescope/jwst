import os

import pytest
from astropy.io.fits.diff import FITSDiff
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_detector1 pipeline on MIRI imaging data."""
    rtdata = rtdata_module
    rtdata.get_data("miri/image/jw00001001001_01101_00001_MIRIMAGE_uncal.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_detector1.cfg", rtdata.input, 
            "--steps.dq_init.save_results=True",
            "--steps.lastframe.save_results=True",
            "--steps.firstframe.save_results=True",
            "--steps.saturation.save_results=True",
            "--steps.rscd.save_results=True",
            "--steps.linearity.save_results=True",
            "--steps.dark_current.save_results=True",
            "--steps.refpix.save_results=True",
            "--steps.jump.save_results=True"]

    Step.from_cmdline(args)
    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("output", ['rate', 'rateints', 'linearity', 'rscd',
                                    'dq_init', 'firstframe', 'lastframe',
                                    'saturation', 'dark_current', 'refpix', 'jump'])
def test_miri_image_caldetector1(run_pipeline, request, fitsdiff_default_kwargs, output):
    """
    Regression test of calwebb_detector1 pipeline performed on MIRI data.
    """

    rtdata = run_pipeline
    rtdata.output = "jw00001001001_01101_00001_MIRIMAGE_" + output + ".fits"

    rtdata.get_truth(os.path.join("truth/test_miri_image_caldetector1",
                                  "jw00001001001_01101_00001_MIRIMAGE_" + output + ".fits"))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
