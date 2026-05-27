"""Regression tests for NIRSpec FS trace modeling."""

import os

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth
TRUTH_PATH = "truth/test_nirspec_fs_trace_model"

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


@pytest.fixture(scope="module")
def run_spec2_trace_model_fs(rtdata_module):
    """Run the Spec2Pipeline and create a trace model with oversampling."""
    rtdata = rtdata_module

    # Get the input rate file
    rtdata.get_asn("nirspec/fs/jw01309-o022_20230113t025924_spec2_00001_asn.json")

    # Run the pipeline
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.adaptive_trace_model.skip=false",
        "--steps.adaptive_trace_model.save_results=true",
        "--steps.adaptive_trace_model.oversample=2.0",
        "--steps.resample_spec.pixel_scale_ratio=0.5",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope="module")
def run_spec2_trace_model_bots(rtdata_module):
    """Run the Spec2Pipeline and create a trace model with oversampling."""
    rtdata = rtdata_module

    # Get the input rate file
    rtdata.get_data("nirspec/tso/jw02420001001_04101_00001-first100_nrs1_rateints.fits")

    # Run the pipeline
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.adaptive_trace_model.skip=false",
        "--steps.adaptive_trace_model.save_results=true",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.mark.parametrize(
    "suffix",
    ["cal", "adaptive_trace_model", "s2d", "x1d"],
)
def test_nirspec_fs_spec2_trace_model(
    run_spec2_trace_model_fs, fitsdiff_default_kwargs, suffix, rtdata_module
):
    """Test trace model and oversample on a NIRSpec FS exposure."""

    rtdata = run_spec2_trace_model_fs
    output = f"jw01309022001_04102_00004_nrs2_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "suffix",
    ["calints", "adaptive_trace_model", "x1dints"],
)
def test_nirspec_bots_spec2_trace_model(
    run_spec2_trace_model_bots, fitsdiff_default_kwargs, suffix, rtdata_module
):
    """Test trace model and oversample on a NIRSpec BOTS exposure."""

    rtdata = run_spec2_trace_model_bots
    output = f"jw02420001001_04101_00001-first100_nrs1_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
