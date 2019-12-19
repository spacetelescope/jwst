import os

import pytest
from astropy.io.fits.diff import FITSDiff
from jwst.pipeline import Detector1Pipeline


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_detector1 pipeline on MIRI imaging data."""
    rtdata = rtdata_module
    rtdata.get_data("miri/image/jw00001001001_01101_00001_MIRIMAGE_uncal.fits")

    step = Detector1Pipeline()
    step.lastframe.save_results = True
    step.firstframe.save_results = True
    step.dq_init.save_results = True
    step.saturation.save_results = True
    step.rscd.save_results = True
    step.linearity.save_results = True
    step.dark_current.save_results = True
    step.refpix.save_results = True
    step.jump.save_results = True
    step.save_results = True
    step.run(rtdata.input)
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
