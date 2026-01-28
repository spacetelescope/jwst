"""Regression tests for MIRI MRS trace modeling."""

import os

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth
TRUTH_PATH = "truth/test_miri_mrs_trace_model"

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


@pytest.fixture(scope="module")
def run_spec2_trace_model(rtdata_module):
    """Run the Spec2Pipeline and create a trace model without oversampling."""
    rtdata = rtdata_module

    # Get the input rate file
    rtdata.get_data("miri/mrs/jw01024001001_04101_00001_mirifulong_rate.fits")

    # Run the pipeline
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.adaptive_trace_model.skip=false",
        "--steps.adaptive_trace_model.save_results=true",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope="module")
def run_spec3_oversample(rtdata_module):
    """Run the Spec3Pipeline on association with oversampling."""
    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01024_ifushort_mediumlong_spec3_00001_asn.json")

    args = [
        "calwebb_spec3",
        rtdata.input,
        "--steps.outlier_detection.skip=true",  # skip for speed
        "--steps.pixel_replace.skip=true",  # skip for speed
        "--steps.adaptive_trace_model.skip=false",
        "--steps.adaptive_trace_model.save_results=true",
        "--steps.adaptive_trace_model.oversample=2.0",
        "--steps.cube_build.save_results=true",
        "--steps.extract_1d.save_results=true",
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.mark.parametrize(
    "suffix",
    ["cal", "adaptive_trace_model", "s3d", "x1d"],
)
def test_miri_mrs_spec2_trace_model(
    run_spec2_trace_model, fitsdiff_default_kwargs, suffix, rtdata_module
):
    """Regression test of the calwebb_spec2 pipeline on a MIRI MRS long-wave exposure"""

    rtdata = run_spec2_trace_model
    output = f"jw01024001001_04101_00001_mirifulong_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "output",
    [
        "jw01024009001_02101_00001_mirifushort_c1000_adaptive_trace_model.fits",
        "jw01024009001_02101_00002_mirifushort_c1000_adaptive_trace_model.fits",
        "jw01024009001_02101_00003_mirifushort_c1000_adaptive_trace_model.fits",
        "jw01024009001_02101_00004_mirifushort_c1000_adaptive_trace_model.fits",
        "jw01024013001_02101_00001_mirifushort_c1000_adaptive_trace_model.fits",
        "jw01024013001_02101_00002_mirifushort_c1000_adaptive_trace_model.fits",
        "jw01024013001_02101_00003_mirifushort_c1000_adaptive_trace_model.fits",
        "jw01024013001_02101_00004_mirifushort_c1000_adaptive_trace_model.fits",
        "jw01024-c1000_t002_miri_ch1-long_x1d.fits",
        "jw01024-c1000_t002_miri_ch1-medium_x1d.fits",
        "jw01024-c1000_t002_miri_ch2-long_x1d.fits",
        "jw01024-c1000_t002_miri_ch2-medium_x1d.fits",
        "jw01024-c1000_t002_miri_ch1-long_s3d.fits",
        "jw01024-c1000_t002_miri_ch1-medium_s3d.fits",
        "jw01024-c1000_t002_miri_ch2-long_s3d.fits",
        "jw01024-c1000_t002_miri_ch2-medium_s3d.fits",
    ],
)
def test_miri_mrs_spec3_oversample(run_spec3_oversample, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""

    rtdata = run_spec3_oversample
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
