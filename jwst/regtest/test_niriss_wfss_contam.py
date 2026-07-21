"""Regression test for WFSS contam correction of NIRISS data"""

import pytest

from jwst.regtest import regtestdata as rt
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_wfss_contam(rtdata_module, resource_tracker):
    """Run the wfss_contam step"""
    rtdata = rtdata_module

    # Get input data files
    rtdata.get_data("niriss/wfss/jw02738-o001_t001_niriss_clear-f200w_i2d.fits")
    rtdata.get_data("niriss/wfss/jw02738-o001_t001_niriss_clear-f200w_segm.fits")
    rtdata.get_data("niriss/wfss/jw02738-o001_t001_niriss_clear-f200w_cat.ecsv")
    rtdata.get_data("niriss/wfss/jw02738001001_04201_00004_nis_rate.fits")
    rtdata.get_data("niriss/wfss/jw02738-o001_20260704t072252_spec2_00013_asn.json")

    # Run the step
    step_params = {
        "step": "jwst.pipeline.Spec2Pipeline",
        "args": [
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


def test_log_tracked_resources_nis_wfss_contam(log_tracked_resources, run_wfss_contam):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["simul", "simul_slits", "contam", "wfss_contam"])
def test_niriss_wfss_contam(run_wfss_contam, fitsdiff_default_kwargs, suffix):
    """Regression test for wfss_contam applied to NIRISS WFSS data"""

    rtdata = run_wfss_contam
    output = f"jw02738001001_04201_00004_nis_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_niriss_wfss_contam/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
