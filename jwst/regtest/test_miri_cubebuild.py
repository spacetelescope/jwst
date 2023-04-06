"""Test CubeBuildStep on MIRI MRS"""
import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope='module')
def run_cube_build_single_output(rtdata_module):
    """Run cube_build on multiple inputs but single output"""
    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/jw01024-o001_20220501t155404_spec3_002_asn.json')

    args = [
        'jwst.cube_build.CubeBuildStep',
        rtdata.input,
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    'output',
    ['jw01024-o001_t002_miri_ch1-short_s3d.fits', 'jw01024-o001_t002_miri_ch2-short_s3d.fits']
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
def test_cube_build_miri_ifualign(rtdata, fitsdiff_default_kwargs):
    """Run cube_build on single file using coord system = ifu_align"""
    input_file = 'jw01024001001_04101_00001_mirifushort_cal.fits'
    rtdata.get_data(f"miri/mrs/{input_file}")

    args = [
        'jwst.cube_build.CubeBuildStep',
        input_file,
        '--coord_system=ifualign'
    ]
    Step.from_cmdline(args)

    output = input_file.replace('cal', 'ch1-short_s3d')
    rtdata.output = output

    rtdata.get_truth(f"truth/test_miri_cubebuild/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
