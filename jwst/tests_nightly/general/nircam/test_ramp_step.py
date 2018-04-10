import os
from astropy.io import fits as pf
from jwst.ramp_fitting.ramp_fit_step import RampFitStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_ramp_fit_nircam():
    """

    Regression test of ramp_fit step performed on NIRCam data.

    """
    output_file_base, output_files = add_suffix('rampfit_output.fits', 'rampfit', list(range(1)))

    try:
        os.remove(output_files[0])
        os.remove("rampfit_opt_out.fits")
    except:
        pass


    RampFitStep.call(BIGDATA+'/nircam/test_ramp_fit/jw00017001001_01101_00001_NRCA1_jump.fits',
                     output_file=output_file_base,
                     save_opt=True,
                     opt_name='rampfit_opt_out.fits'
                     )

    # compare primary output
    n_priout = output_files[0]
    h = pf.open( n_priout )
    n_priref = BIGDATA+'/nircam/test_ramp_fit/jw00017001001_01101_00001_NRCA1_ramp_fit.fits'
    href = pf.open( n_priref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the standard output file - a:', n_priout)
    print (' ... and the reference file - b:', n_priref)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)


    # compare optional output
    n_optout = 'rampfit_opt_out_fitopt.fits'
    h = pf.open( n_optout)
    n_optref = BIGDATA+'/nircam/test_ramp_fit/rampfit_opt_out.fits'
    href = pf.open( n_optref )
    newh = pf.HDUList([h['primary'],h['slope'],h['sigslope'],h['yint'],h['sigyint'],h['pedestal'],h['weights'],h['crmag']])
    newhref = pf.HDUList([href['primary'],href['slope'],href['sigslope'],href['yint'],href['sigyint'],href['pedestal'],href['weights'],href['crmag']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the optional output file - a:', n_optout)
    print (' ... and the reference file - b:', n_optref)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)
