import pytest
from astropy.io import fits
from jwst.pipeline.calwebb_dark import DarkPipeline

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_nirspec_dark_pipeline(_bigdata):
    """
    Regression test of calwebb_dark pipeline performed on NIRSpec raw data.
    """
    pipe = DarkPipeline()
    pipe.suffix = 'dark'
    pipe.ipc.skip = True
    pipe.refpix.odd_even_columns = True
    pipe.refpix.use_side_ref_pixels = True
    pipe.refpix.side_smoothing_length = 11
    pipe.refpix.side_gain = 1.0
    pipe.refpix.odd_even_rows = True
    pipe.output_file = 'jw84500013001_02103_00003_NRS1_uncal.fits'

    pipe.run(_bigdata+'/pipelines/jw84500013001_02103_00003_NRS1_uncal.fits')

    h = fits.open('jw84500013001_02103_00003_NRS1_dark.fits')
    href = fits.open(_bigdata+'/pipelines/jw84500013001_02103_00003_NRS1_dark_ref.fits')
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = fits.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)
    assert result.identical, result.report()
