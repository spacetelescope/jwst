import pytest
from astropy.io import fits
from jwst.cube_build.cube_build_step import CubeBuildStep
from jwst import datamodels

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_cubebuild_miri(_bigdata):
    """
    Regression test of cube_build performed on MIRI MRS data.
    """
    input_model = datamodels.IFUImageModel(_bigdata+'/miri/test_cube_build/jw10001001001_01101_00001_mirifushort_cal.fits')
    CubeBuildStep.call(input_model, output_type='multi', save_results=True)

    h = fits.open('jw10001001001_01101_00001_mirifushort_s3d.fits')
    href = fits.open(_bigdata+'/miri/test_cube_build/jw10001001001_01101_00001_mirifushort_s3d_ref.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['wmap']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['wmap']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001
    )
    assert result.identical, result.report()
