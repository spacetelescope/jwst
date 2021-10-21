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
    rtdata.get_truth(f'truth/test_miri_cubebuild/{output}')

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_cube_build_miri_internal_cal(rtdata, fitsdiff_default_kwargs):
    """Run cube_build on single file using coord system = internal_cal"""
    input_file = 'det_image_seq2_MIRIFUSHORT_12SHORTexp1_cal.fits'
    rtdata.get_data(f"miri/mrs/{input_file}")

    args = [
        'jwst.cube_build.CubeBuildStep',
        input_file,
        '--save_results=true',
        '--coord_system=internal_cal'
    ]
    Step.from_cmdline(args)

    output = input_file.replace('cal', 'ch1-short_internal_s3d')
    rtdata.output = output

    rtdata.get_truth(f"truth/test_miri_cubebuild/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
