"""Regression tests for Moving Target MIRI MRS mode"""
import os
import pytest
from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step

# Define artifactory source and truth

TRUTH_PATH = 'truth/test_miri_mrs_spec3_moving_target/'


@pytest.fixture(scope='module')
def run_spec3_moving_target(jail, rtdata_module):
    """Run the Spec3Pipeline dithered flight data """

    # Association has 2 exposures from IFUSHORT

    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/jw01449-o002_ch12_spec3_00001_asn.json')

    args = [
        "calwebb_spec3",
        rtdata.input,
        '--steps.outlier_detection.skip=True',
        '--steps.cube_build.save_results=true',
        '--steps.extract_1d.save_results=true',
    ]

    Step.from_cmdline(args)
    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    [
        'jw01449-o002_t001_miri_ch1-short_s3d.fits',
        'jw01449-o002_t001_miri_ch2-short_s3d.fits',
        'jw01449-o002_t001_miri_ch1-short_x1d.fits',
        'jw01449-o002_t001_miri_ch2-short_x1d.fits'
    ],

)
def test_spec3_moving_target(run_spec3_moving_target, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""

    rtdata = run_spec3_moving_target
    rtdata.output = output

    rtdata.get_truth(os.path.join(TRUTH_PATH, rtdata.output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
