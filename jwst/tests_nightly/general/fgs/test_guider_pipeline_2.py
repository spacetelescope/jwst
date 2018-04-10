import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_guider import GuiderPipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_guider_pipeline2():
    """

    Regression test of calwebb_guider pipeline performed on ACQ-1 data.

    """

    GuiderPipeline.call(BIGDATA+'/fgs/test_guiderpipeline/jw88600073001_gs-acq1_2016022183837_uncal.fits',
                        output_file='jw88600073001_gs-acq1_2016022183837_cal.fits')

    # Compare calibrated ramp product
    n_cr = 'jw88600073001_gs-acq1_2016022183837_cal.fits'
    h = pf.open( n_cr )
    n_ref = BIGDATA+'/fgs/test_guiderpipeline/jw88600073001_gs-acq1_2016022183837_cal_ref.fits'
    href = pf.open( n_ref )
    newh = pf.HDUList([h['primary'],h['sci'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001
    )

    print (' Fitsdiff comparison between the calibrated product file - a:', n_cr )
    print (' ... and the reference file - b:', n_ref)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

