import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_ami3 import Ami3Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_ami_pipeline():
    """

    Regression test of the AMI pipeline performed on NIRISS AMI data.

    """

    Ami3Pipeline.call(BIGDATA+'/niriss/test_ami_pipeline/test_lg1_asn.json',
                      config_file='calwebb_ami3.cfg')
    
    h = pf.open('test_targ_aminorm.fits')
    href = pf.open(BIGDATA+'/niriss/test_ami_pipeline/ami_pipeline_targ_lgnorm.fits') 
    newh = pf.HDUList([h['primary'],h['fit'],h['resid'],h['closure_amp'],
                       h['closure_pha'],h['fringe_amp'],h['fringe_pha'],
                       h['pupil_pha'],h['solns']])
    newhref = pf.HDUList([href['primary'],href['fit'],href['resid'],href['closure_amp'],
                          href['closure_pha'],href['fringe_amp'],href['fringe_pha'],
                          href['pupil_pha'],href['solns']])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

