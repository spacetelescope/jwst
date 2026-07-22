"""Regression test for WFSS contam correction of NIRCam data"""

import pytest

from jwst.regtest import regtestdata as rt
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_spec2_with_contam(rtdata_module, resource_tracker):
    """Run the Spec2Pipeline with the wfss_contam step turned on."""
    rtdata = rtdata_module

    # Get input data files
    rtdata.get_data("nircam/wfss/jw01076-o107_t001_nircam_clear-f322w2_i2d.fits")
    rtdata.get_data("nircam/wfss/jw01076-o107_t001_nircam_clear-f322w2_segm.fits")
    rtdata.get_data("nircam/wfss/jw01076-o107_t001_nircam_clear-f322w2_cat.ecsv")
    rtdata.get_data("nircam/wfss/jw01076107001_02101_00002_nrcalong_rate.fits")
    rtdata.get_data("nircam/wfss/jw01076-o107_20260704t223417_spec2_00006_asn.json")

    # Run the step
    step_params = {
        "step": "jwst.pipeline.Spec2Pipeline",
        "args": [
            "--steps.extract_2d.wfss_nbright=10",
            "--steps.wfss_contam.magnitude_limit=20",
            "--steps.wfss_contam.save_simulated_image=True",
            "--steps.wfss_contam.save_contam_images=True",
            "--steps.wfss_contam.save_results=true",
            "--steps.wfss_contam.skip=False",
            "--steps.wfss_contam.maximum_cores=none",
        ],
    }
    with resource_tracker.track():
        rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


def test_log_tracked_resources_nircam_wfss_contam(log_tracked_resources, run_spec2_with_contam):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["simul", "contam", "wfss_contam"])
def test_nrc_wfss_contam(run_spec2_with_contam, fitsdiff_default_kwargs, suffix):
    """Regression test for wfss_contam applied to NIRCam WFSS data"""
    rtdata = run_spec2_with_contam
    output = f"jw01076107001_02101_00002_nrcalong_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nircam_wfss_contam/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
