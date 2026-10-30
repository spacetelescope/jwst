"""Regression tests for MIRI MRS charge migration correction."""

import os

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth
TRUTH_PATH = "truth/test_miri_mrs_charge_migration"

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


@pytest.fixture(scope="module")
def run_detector1_charge_migration(rtdata_module):
    """Run the Detector1Pipeline with charge migration correction."""
    rtdata = rtdata_module

    # Get the input rate file
    rtdata.get_data("miri/mrs/jw03779002001_03102_00001_mirifulong_uncal.fits")

    # Run the pipeline
    args = [
        "calwebb_detector1",
        rtdata.input,
        "--steps.charge_migration.skip=False",
        "--steps.charge_migration.save_results=True",
        "--steps.charge_migration.correction_threshold=10000",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.mark.parametrize("suffix", ["charge_migration", "rate"])
def test_miri_mrs_det1_charge_migration(
    run_detector1_charge_migration, fitsdiff_default_kwargs, suffix, rtdata_module
):
    """Regression test of the calwebb_spec2 pipeline on a MIRI MRS long-wave exposure"""

    rtdata = run_detector1_charge_migration
    output = f"jw03779002001_03102_00001_mirifulong_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
