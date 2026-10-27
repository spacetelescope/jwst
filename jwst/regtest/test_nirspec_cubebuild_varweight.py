"""Test CubeBuildStep on NIRSpec data using variance weighting"""

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


@pytest.fixture(scope="module")
def run_cubebuild_nirspec_readnoise(rtdata_module, resource_tracker):
    """Run cube_build on nirspec data using drizzle and readnoise weighting"""
    input_file = "jw01794-o001_20261005t183207_spec3_00001_1set.asn.json"
    # Test has F100LP and G140M
    rtdata = rtdata_module

    rtdata.get_asn(f"nirspec/ifu/{input_file}")

    args = ["jwst.cube_build.CubeBuildStep", rtdata.input, "--variance_weighting=readnoise"]

    with resource_tracker.track():
        Step.from_cmdline(args)
    return rtdata


def test_log_tracked_resources(log_tracked_resources, run_cubebuild_nirspec_readnoise):
    log_tracked_resources()


@pytest.mark.parametrize(
    "output",
    ["jw01794-o001_t001_nirspec_g140m-f100lp_s3d.fits"],
)
def test_cubebuild_nirspec_readnoise(
    run_cubebuild_nirspec_readnoise, fitsdiff_default_kwargs, output
):
    rtdata = run_cubebuild_nirspec_readnoise
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_cubebuild_varweight/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
