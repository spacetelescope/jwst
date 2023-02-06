"""Regression tests for MIRI MRS modes"""
import os
import pytest
from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth
TRUTH_PATH = 'truth/test_miri_mrs'


@pytest.fixture(scope='module')
def run_spec3_single(jail, rtdata_module):
    """Run the Spec3Pipeline on single cal from Spec2Pipeline"""

    # Running Spec2 on MIRI data creates a single cube with 2 channels
    # Running Spec3 on MIRI data creates 2 channel cubes
    # This test will perform a regression tests on the 2 channel cubes

    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/jw01024-o001_single_spec3_00005_asn.json')

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.cube_build.save_results=true',
        '--steps.extract_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata


@pytest.fixture(scope='module')
def run_spec3_multi(jail, rtdata_module):
    """Run the Spec3Pipeline dithered flight data """

    # Use an association that has at least 2 bands and both IFULONG & IFUSHORT
    # so we can create multi band channel cubes from the two detectors.
    # We will run tests on a cubes from IFUSHORT and  from  IFULONG

    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/jw01024-c10000_spec3_00001_asn.json')

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.master_background.save_results=true',
        '--steps.outlier_detection.save_results=true',
        '--steps.cube_build.save_results=true',
        '--steps.extract_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw01024001001_04101_00001_mirifulong_ch3-short_s3d.fits',
        'jw01024001001_04101_00001_mirifulong_ch3-short_x1d.fits',
        'jw01024001001_04101_00001_mirifulong_ch4-short_s3d.fits',
        'jw01024001001_04101_00001_mirifulong_ch4-short_x1d.fits',
    ],
)
def test_spec3_single(run_spec3_single, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rtdata = run_spec3_single
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, rtdata.output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw01024001001_04101_00001_mirifushort_c1000_crf.fits',
        'jw01024001001_04101_00002_mirifushort_c1000_crf.fits',
        'jw01024001001_04101_00001_mirifulong_c1000_crf.fits',
        'jw01024001001_04101_00002_mirifulong_c1000_crf.fits',
        'jw01024-c1000_miri_ch1-shortmedium_s3d.fits',
        'jw01024-c1000_miri_ch2-shortmedium_s3d.fits',
        'jw01024-c1000_miri_ch3-shortmedium_s3d.fits',
        'jw01024-c1000_miri_ch4-shortmedium_s3d.fits',
        'jw01024-c1000_miri_ch1-shortmedium_x1d.fits',
        'jw01024-c1000_miri_ch2-shortmedium_x1d.fits',
        'jw01024-c1000_miri_ch3-shortmedium_x1d.fits',
        'jw01024-c1000_miri_ch4-shortmedium_x1d.fits'
    ],
)
def test_spec3_multi(run_spec3_multi, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""

    rtdata = run_spec3_multi
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, rtdata.output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
