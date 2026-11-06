"""Test CubeBuildStep on MIRI data using variance weighting"""

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


@pytest.fixture(scope="module")
def run_cubebuild_miri_readnoise(rtdata_module, resource_tracker):
    """Run cube_build on miri Channel 2 Medium  data using drizzle and readnoise weighting"""
    input_file = "jw01523-o004_20260801t122412_spec3_00001_ch12_asn.json"
    # Test has F100LP and G140M
    rtdata = rtdata_module

    rtdata.get_asn(f"miri/mrs/{input_file}")

    args = [
        "jwst.cube_build.CubeBuildStep",
        rtdata.input,
        "--variance_weighting=readnoise",
        "--channel=2",
    ]

    with resource_tracker.track():
        Step.from_cmdline(args)
    return rtdata


def test_log_tracked_resources(log_tracked_resources, run_cubebuild_miri_readnoise):
    log_tracked_resources()


@pytest.mark.parametrize(
    "output",
    ["jw01523-o004_t009_miri_ch2-medium_s3d.fits"],
)
def test_cubebuild_miri_readnoise(run_cubebuild_miri_readnoise, fitsdiff_default_kwargs, output):
    rtdata = run_cubebuild_miri_readnoise
    rtdata.output = output

    rtdata.get_truth(f"truth/test_miri_cubebuild_varweight/{output}")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
