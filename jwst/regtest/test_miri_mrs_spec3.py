"""Regression tests for MIRI MRS modes"""

import os

import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

# Define artifactory source and truth
TRUTH_PATH = "truth/test_miri_mrs"

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


@pytest.fixture(scope="module")
def run_spec3_ifushort(rtdata_module, resource_tracker):
    """Run the Spec3Pipeline on association with 2 bands on IFUSHORT"""

    # Test has bands medium and long for IFUSHORT

    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01024_ifushort_mediumlong_spec3_00001_asn.json")

    args = [
        "calwebb_spec3",
        rtdata.input,
        "--steps.outlier_detection.save_results=true",
        "--steps.cube_build.save_results=true",
        "--steps.extract_1d.save_results=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope="module")
def run_spec3_ifulong(rtdata_module, resource_tracker):
    """Run the Spec3Pipeline dithered flight data"""

    # Test has bands medium and long for IFULONG

    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01024_ifulong_mediumlong_spec3_00001_asn.json")

    args = [
        "calwebb_spec3",
        rtdata.input,
        "--steps.cube_build.save_results=true",
        "--steps.extract_1d.save_results=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope="module")
def run_spec3_ifushort_emsm(rtdata_module):
    """Run the Spec3Pipeline (cube_build using weighting emsm) on association with 2 bands on IFUSHORT"""

    # Test has bands medium and long for IFUSHORT

    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01024_ifushort_mediumlong_spec3_00001_asn.json")

    args = [
        "calwebb_spec3",
        rtdata.input,
        "--steps.cube_build.save_results=true",
        "--steps.cube_build.weighting=emsm",
        '--steps.cube_build.output_file="miri_mrs_emsm"',
        "--steps.extract_1d.save_results=true",
    ]
    with pytest.warns(Warning, match="Sources were found, but none pass"):
        Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope="module")
def run_spec3_ifushort_extract1d(rtdata_module):
    """Run the Spec3Pipeline on association with 2 bands on IFUSHORT"""

    # Test has bands medium and long for IFUSHORT

    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01024_ifushort_mediumlong_spec3_extract1d_00001_asn.json")

    args = [
        "calwebb_spec3",
        rtdata.input,
        "--steps.outlier_detection.save_results=true",
        "--steps.cube_build.save_results=true",
        '--steps.extract_1d.ifu_set_srctype="POINT"',
        "--steps.extract_1d.ifu_rscale=3.0",
        "--steps.extract_1d.ifu_rfcorr=true",
        "--steps.extract_1d.save_results=true",
    ]
    Step.from_cmdline(args)
    return rtdata


def test_log_tracked_resources_spec3short(log_tracked_resources, run_spec3_ifushort):
    log_tracked_resources()


def test_log_tracked_resources_spec3long(log_tracked_resources, run_spec3_ifulong):
    log_tracked_resources()


@pytest.mark.parametrize(
    "output",
    [
        "jw01024-c1000_t002_miri_ch3-long_x1d.fits",
        "jw01024-c1000_t002_miri_ch3-medium_x1d.fits",
        "jw01024-c1000_t002_miri_ch4-long_x1d.fits",
        "jw01024-c1000_t002_miri_ch4-medium_x1d.fits",
        "jw01024-c1000_t002_miri_ch3-long_s3d.fits",
        "jw01024-c1000_t002_miri_ch3-medium_s3d.fits",
        "jw01024-c1000_t002_miri_ch4-long_s3d.fits",
        "jw01024-c1000_t002_miri_ch4-medium_s3d.fits",
    ],
)
def test_spec3_ifulong(run_spec3_ifulong, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rtdata = run_spec3_ifulong
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "output",
    [
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
def test_spec3_ifushort(run_spec3_ifushort, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""

    rtdata = run_spec3_ifushort
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "output",
    [
        "miri_mrs_emsm_ch1-long_x1d.fits",
        "miri_mrs_emsm_ch1-medium_x1d.fits",
        "miri_mrs_emsm_ch2-long_x1d.fits",
        "miri_mrs_emsm_ch2-medium_x1d.fits",
        "miri_mrs_emsm_ch1-long_s3d.fits",
        "miri_mrs_emsm_ch1-medium_s3d.fits",
        "miri_mrs_emsm_ch2-long_s3d.fits",
        "miri_mrs_emsm_ch2-medium_s3d.fits",
    ],
)
def test_spec3_ifushort_emsm(run_spec3_ifushort_emsm, fitsdiff_default_kwargs, output):
    """Regression test using weighting = 'emsm'"""

    rtdata = run_spec3_ifushort_emsm
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "output",
    [
        "jw01024-c1000_t002_extract1dtest_miri_ch1-long_x1d.fits",
        "jw01024-c1000_t002_extract1dtest_miri_ch1-medium_x1d.fits",
        "jw01024-c1000_t002_extract1dtest_miri_ch2-long_x1d.fits",
        "jw01024-c1000_t002_extract1dtest_miri_ch2-medium_x1d.fits",
    ],
)
def test_spec3_ifushort_extract1d(run_spec3_ifushort_extract1d, fitsdiff_default_kwargs, output):
    """Regression test for extract_1d using ifu_set_srctype=POINT, ifu_rscale=3.0,
    ifu_rfcorr=true"""

    rtdata = run_spec3_ifushort_extract1d
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
