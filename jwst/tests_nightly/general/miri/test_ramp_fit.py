import os
from astropy.io import fits as pf
from jwst.ramp_fitting.ramp_fit_step import RampFitStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_ramp_fit_miri1():
    """

    Regression test of ramp_fit step performed on MIRI data.

    """
    output_file_base, output_files = add_suffix('rampfit1_output.fits', 'rampfit', list(range(2)))

    try:
        for output_file in output_files:
            os.remove(output_file)
        os.remove("rampfit1_opt_out_fitopt.fits")
    except:
        pass

    RampFitStep.call(BIGDATA+'/miri/test_ramp_fit/jw00001001001_01101_00001_MIRIMAGE_jump.fits',
                      config_file='ramp_fit.cfg',
                      save_opt=True,
                      opt_name='rampfit1_opt_out.fits',
                      output_file=output_file_base
    )

    # compare primary output
    n_priout = output_files[0]
    h = pf.open( n_priout )
    n_priref = BIGDATA+'/miri/test_ramp_fit/jw00001001001_01101_00001_MIRIMAGE_ramp_fit.fits'
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

    # compare integration-specific output
    n_intout = output_files[1]
    h = pf.open( n_intout )
    n_intref = BIGDATA+'/miri/test_ramp_fit/jw00001001001_01101_00001_MIRIMAGE_int.fits'
    href = pf.open( n_intref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )

    print (' Fitsdiff comparison between the integration-specfic output file - a:', n_intout)
    print (' ... and the reference file - b:', n_intref)

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    # compare optional output
    n_optout = 'rampfit1_opt_out_fitopt.fits'
    h = pf.open( n_optout )
    n_optref = BIGDATA+'/miri/test_ramp_fit/jw00001001001_01101_00001_MIRIMAGE_opt.fits'
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
