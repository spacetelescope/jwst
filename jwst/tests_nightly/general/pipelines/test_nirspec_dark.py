import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_dark import DarkPipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_nirspec_dark_pipeline():
    """

    Regression test of calwebb_dark pipeline performed on NIRSpec raw data.

    """
    step = DarkPipeline()
    step.suffix = "dark"
    step.refpix.odd_even_columns = True
    step.refpix.use_side_ref_pixels = True
    step.refpix.side_smoothing_length=11
    step.refpix.side_gain=1.0

    step.run(BIGDATA+'/pipelines/jw84500013001_02103_00003_NRS1_uncal.fits',
             output_file='jw84500013001_02103_00003_NRS1_dark.fits')

    h = pf.open('jw84500013001_02103_00003_NRS1_dark.fits')
    href = pf.open(BIGDATA+'/pipelines/jw84500013001_02103_00003_NRS1_dark_ref.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

