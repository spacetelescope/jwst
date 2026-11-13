"""Regression tests for MIRI MRS straylight step"""

import pytest

from jwst.assign_wcs import AssignWcsStep
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.straylight import StraylightStep

# Define artifactory source and truth file paths
INPUT_PATH = "miri/mrs"
TRUTH_PATH = "truth/test_miri_mrs_straylight"


@pytest.fixture(scope="module")
def run_straylight_clean_showers(rtdata_module):
    """Run straylight step turning on clean shower and save_shower model."""
    rtdata = rtdata_module

    # Get the MIRI MRS rate file
    filename = "jw01264014001_02101_00001_mirifushort_rate.fits"
    rtdata.get_data(INPUT_PATH + "/" + filename)

    with AssignWcsStep.call(rtdata.input) as wcs_model:
        result = StraylightStep.call(
            wcs_model,
            clean_showers=True,
            save_shower_model=True,
            save_results=True,
        )
        result.close()

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "suffix",
    ["straylightstep", "shower_model"],
)
def test_miri_mrs_straylight_clean_showers(run_straylight_clean_showers, fitsdiff_default_kwargs, suffix):
    """Test running straylight with clean shower= True on an rate file with significant showers."""

    rtdata = run_straylight_clean_showers
    output = "jw01264014001_02101_00001_mirifushort_" + suffix + ".fits"
    rtdata.output = output
    # Get the truth file
    rtdata.get_truth(f"{TRUTH_PATH}/{output}")

    # Compare the results for straylight output with clean_showers turned on.
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.fixture(scope="module")
def run_straylight_xartifact(rtdata_module):
    """Run straylight step with defaults."""
    rtdata = rtdata_module

    # Get the MIRI MRS rate file
    filename = "jw01524017001_04101_00002_mirifushort_rate.fits"
    rtdata.get_data(INPUT_PATH + "/" + filename)

    with AssignWcsStep.call(rtdata.input) as wcs_model:
        result = StraylightStep.call(
            wcs_model,
            save_results=True,
        )
        result.close()

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "suffix",
    ["straylightstep"],
)
def test_miri_mrs_straylight_xartifact(run_straylight_xartifact, fitsdiff_default_kwargs, suffix):
    """Test running straylight on a rate file with significant cross-artifact."""

    rtdata = run_straylight_xartifact
    output = "jw01524017001_04101_00002_mirifushort_" + suffix + ".fits"
    rtdata.output = output
    # Get the truth file
    rtdata.get_truth(f"{TRUTH_PATH}/{output}")

    # Compare the results for straylight output with clean_showers turned on.
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
