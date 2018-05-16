import os
import pytest
from astropy.io import fits as pf
from jwst.ramp_fitting.ramp_fit_step import RampFitStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_ramp_fit_niriss(_bigdata):
    """

    Regression test of ramp_fit step performed on NIRISS data.

    """
    suffix = 'rampfit'
    output_file_base, output_files = add_suffix('rampfit_output.fits', suffix, list(range(1)))

    try:
        for output_file in output_files:
            os.remove(output_file)
        os.remove("rampfit_opt_out.fits")
    except:
        pass


    RampFitStep.call(_bigdata+'/niriss/test_ramp_fit/jw00034001001_01101_00001_NIRISS_jump.fits',
                      save_opt=True, output_file=output_file_base, suffix=suffix,
                      opt_name='rampfit_opt_out.fits'
    )
    # primary output
    n_priout = output_files[0]
    h = pf.open( n_priout )
    n_priref = _bigdata+'/niriss/test_ramp_fit/jw00034001001_01101_00001_NIRISS_ramp_fit.fits'
    href = pf.open( n_priref)
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()


    # optional output
    n_optout = 'rampfit_opt_out_fitopt.fits'
    h = pf.open( n_optout )
    n_optref = _bigdata+'/niriss/test_ramp_fit/jw00034001001_01101_00001_NIRISS_uncal_opt.fits'
    href = pf.open( n_optref )
    newh = pf.HDUList([h['primary'],h['slope'],h['sigslope'],h['yint'],h['sigyint'],h['pedestal'],h['weights'],h['crmag']])
    newhref = pf.HDUList([href['primary'],href['slope'],href['sigslope'],href['yint'],href['sigyint'],href['pedestal'],href['weights'],href['crmag']])

    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
