"""Test CubeBuildStep on MIRI MRS"""
import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope='module')
def run_cube_build_single_output(rtdata_module):
    """Run cube_build on multiple inputs but single output"""
    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/two_spec3_asn.json')

    args = [
        'jwst.cube_build.CubeBuildStep',
        rtdata.input,
        '--save_results=true'
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.slow
@pytest.mark.parametrize(
    'output',
    ['outtest_ch1-short_s3d.fits', 'outtest_ch2-short_s3d.fits']
)
def test_cube_build_single_output(run_cube_build_single_output, output, fitsdiff_default_kwargs):
    """Test just running cube build and ensure that output happens"""
    rtdata = run_cube_build_single_output
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth('truth/test_miri_cubebuild/' +  output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
