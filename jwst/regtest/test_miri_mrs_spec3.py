"""Regression tests for MIRI MRS modes"""
import os
import pytest
from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth
TRUTH_PATH = 'truth/test_miri_mrs'


@pytest.fixture(scope='module')
def run_spec3_ifushort(jail, rtdata_module):
    """Run the Spec3Pipeline on association with 2 bands on IFUSHORT"""

    # Test has bands medium and long for IFUSHORT

    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/jw01024_ifushort_mediumlong_spec3_00001_asn.json')

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.outlier_detection.save_results=true',
        '--steps.cube_build.save_results=true',
        '--steps.extract_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope='module')
def run_spec3_ifulong(jail, rtdata_module):
    """Run the Spec3Pipeline dithered flight data """

    # Test has bands medium and long for IFULONG

    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/jw01024_ifulong_mediumlong_spec3_00001_asn.json')

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.cube_build.save_results=true',
        '--steps.extract_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope='module')
def run_spec3_ifushort_emsm(jail, rtdata_module):
    """Run the Spec3Pipeline (cube_build using weighting emsm) on association with 2 bands on IFUSHORT"""

    # Test has bands medium and long for IFUSHORT

    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/jw01024_ifushort_mediumlong_spec3_00001_asn.json')

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.cube_build.save_results=true',
        '--steps.cube_build.weighting=emsm',
        '--steps.cube_build.output_file="miri_mrs_emsm"',
        '--steps.extract_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope='module')
def run_spec3_ifushort_extract1d(jail, rtdata_module):
    """Run the Spec3Pipeline on association with 2 bands on IFUSHORT"""

    # Test has bands medium and long for IFUSHORT

    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/jw01024_ifushort_mediumlong_spec3_extract1d_00001_asn.json')

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.outlier_detection.save_results=true',
        '--steps.cube_build.save_results=true',
        '--steps.extract_1d.ifu_set_srctype="POINT"',
        '--steps.extract_1d.ifu_rscale=3.0',
        '--steps.extract_1d.ifu_rfcorr=true',
        '--steps.extract_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw01024-c1000_t002_miri_ch3-mediumlong_x1d.fits',
        'jw01024-c1000_t002_miri_ch4-mediumlong_x1d.fits',
        'jw01024-c1000_t002_miri_ch3-mediumlong_s3d.fits',
        'jw01024-c1000_t002_miri_ch4-mediumlong_s3d.fits',
        ]
)
def test_spec3_ifulong(run_spec3_ifulong, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rtdata = run_spec3_ifulong
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, rtdata.output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw01024-c1000_t002_miri_ch1-mediumlong_x1d.fits',
        'jw01024-c1000_t002_miri_ch2-mediumlong_x1d.fits',
        'jw01024-c1000_t002_miri_ch1-mediumlong_s3d.fits',
        'jw01024-c1000_t002_miri_ch2-mediumlong_s3d.fits'
    ],
)
def test_spec3_ifushort(run_spec3_ifushort, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""

    rtdata = run_spec3_ifushort
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, rtdata.output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'miri_mrs_emsm_ch1-mediumlong_x1d.fits',
        'miri_mrs_emsm_ch2-mediumlong_x1d.fits',
        'miri_mrs_emsm_ch1-mediumlong_s3d.fits',
        'miri_mrs_emsm_ch2-mediumlong_s3d.fits'
    ],
)
def test_spec3_ifushort_emsm(run_spec3_ifushort_emsm, fitsdiff_default_kwargs, output):
    """Regression test using weighting = 'emsm' """

    rtdata = run_spec3_ifushort_emsm
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, rtdata.output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw01024-c1000_t002_extract1dtest_miri_ch1-mediumlong_x1d.fits',
        'jw01024-c1000_t002_extract1dtest_miri_ch2-mediumlong_x1d.fits'
    ],
)
def test_spec3_ifushort_extract1d(run_spec3_ifushort_extract1d, fitsdiff_default_kwargs, output):
    """Regression test for extract_1d using ifu_set_srctype=POINT, ifu_rscale=3.0,
    ifu_rfcorr=true"""

    rtdata = run_spec3_ifushort_extract1d
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, rtdata.output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
