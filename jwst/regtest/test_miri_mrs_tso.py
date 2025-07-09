"""Regression test for MIRI MRS TSO mode"""

import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

# Define artifactory source and truth
INPUT_PATH = "miri/mrs"
TRUTH_PATH = "truth/test_miri_mrs_tso"

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_spec2(rtdata_module, resource_tracker):
    """Run the Spec2Pipeline on a single exposure"""
    rtdata = rtdata_module

    # Setup the inputs
    file_name = "jw01556001001_04102_00001-seg001_mirifushort_rateints.fits"
    rtdata.get_data(INPUT_PATH + "/" + file_name)

    # Run the pipeline
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.flat_field.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.fringe.save_results=true",
        "--steps.photom.save_results=true",
        "--steps.photom.mrs_time_correction=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_spec2(log_tracked_resources, run_spec2):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix", ["assign_wcs", "calints", "flat_field", "fringe", "photom", "srctype"]
)
def test_spec2(rtdata_module, run_spec2, fitsdiff_default_kwargs, suffix):
    """Test ensuring the calwebb_tso-spec2 is operating appropriately for MIRI MRS TSO data"""
    rtdata = rtdata_module
    output = f"jw01556001001_04102_00001-seg001_mirifushort_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"{TRUTH_PATH}/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
