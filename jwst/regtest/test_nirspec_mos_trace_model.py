"""Regression tests for NIRSpec MOS trace modeling."""

import os

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth
TRUTH_PATH = "truth/test_nirspec_mos_trace_model"

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


@pytest.fixture(scope="module")
def run_spec3_trace_model_mos(rtdata_module):
    """Run the Spec3Pipeline and create a trace model with oversampling."""
    rtdata = rtdata_module

    # Get the MSA metadata file referenced in the input exposure
    rtdata.get_data("nirspec/mos/jw01345066001_01_msa.fits")

    # Get the input ASN file and exposures
    rtdata.get_asn("nirspec/mos/jw01345-o066_20230831t181155_spec3_00002_asn.json")

    # Run the pipeline
    args = [
        "calwebb_spec3",
        rtdata.input,
        "--steps.adaptive_trace_model.skip=false",
        "--steps.adaptive_trace_model.save_results=true",
        "--steps.adaptive_trace_model.oversample=2.0",
        "--steps.resample_spec.pixel_scale_ratio=0.5",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.mark.parametrize(
    "suffix",
    ["adaptive_trace_model", "s2d", "x1d"],
)
def test_nirspec_mos_spec3_trace_model(
    run_spec3_trace_model_mos, fitsdiff_default_kwargs, suffix, rtdata_module
):
    """Test trace model and oversample on a NIRSpec MOS exposure."""
    # Check only one of the sources: the rest are background, faint, or virtual
    # and do not have non-trivial trace models.
    rtdata = run_spec3_trace_model_mos
    output = f"jw01345-o066_s000007380_nirspec_f170lp-g235m_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
