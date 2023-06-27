"""Test CubeBuildStep NIRSPEC internal_cal cubes"""
import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.mark.bigdata
def test_cube_build_nirspec_internal_cal(rtdata, fitsdiff_default_kwargs):
    """Run cube_build on single file using coord system = internal_cal"""
    input_file = 'jw01249005001_03101_00004_nrs1_cal.fits'
    rtdata.get_data(f"nirspec/ifu/{input_file}")

    args = [
        'jwst.cube_build.CubeBuildStep',
        input_file,
        '--save_results=true',
        '--coord_system=internal_cal'
    ]
    Step.from_cmdline(args)

    output = input_file.replace('cal', 'g395h-f290lp_internal_s3d')
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_ifu_internal/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
