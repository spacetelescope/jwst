"""Regression tests for NIRSpec trace modeling."""

import os

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth
TRUTH_PATH = "truth/test_nirspec_ifu_trace_model"

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


@pytest.fixture(scope="module")
def run_spec2_trace_model(rtdata_module):
    """Run the Spec2Pipeline and create a trace model with oversampling."""
    rtdata = rtdata_module

    # Get the input rate file
    rtdata.get_data("nirspec/ifu/jw01251004001_03107_00001_nrs1_rate.fits")

    # Run the pipeline
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.adaptive_trace_model.skip=false",
        "--steps.adaptive_trace_model.save_results=true",
        "--steps.adaptive_trace_model.oversample=2.0",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.mark.parametrize(
    "suffix",
    ["cal", "adaptive_trace_model", "s3d", "x1d"],
)
def test_nirspec_ifu_spec2_trace_model(
    run_spec2_trace_model, fitsdiff_default_kwargs, suffix, rtdata_module
):
    """Test trace model and oversample on a NIRSpec IFU exposure."""

    rtdata = run_spec2_trace_model
    output = f"jw01251004001_03107_00001_nrs1_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
