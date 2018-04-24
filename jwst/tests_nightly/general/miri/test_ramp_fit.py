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

def test_ramp_fit_miri1(_bigdata):
    """

    Regression test of ramp_fit step performed on MIRI data.

    """
    suffix = 'rampfit'
    output_file_base, output_files = add_suffix('rampfit1_output.fits', suffix, list(range(2)))

    try:
        for output_file in output_files:
            os.remove(output_file)
        os.remove("rampfit1_opt_out_fitopt.fits")
    except:
        pass

    RampFitStep.call(_bigdata+'/miri/test_ramp_fit/jw00001001001_01101_00001_MIRIMAGE_jump.fits',
                     save_opt=True,
                     opt_name='rampfit1_opt_out.fits',
                     output_file=output_file_base, suffix=suffix
                     )

    # compare primary output
    n_priout = output_files[0]
    h = pf.open( n_priout )
    n_priref = _bigdata+'/miri/test_ramp_fit/jw00001001001_01101_00001_MIRIMAGE_ramp_fit.fits'
    href = pf.open( n_priref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # compare integration-specific output
    n_intout = output_files[1]
    h = pf.open( n_intout )
    n_intref = _bigdata+'/miri/test_ramp_fit/jw00001001001_01101_00001_MIRIMAGE_int.fits'
    href = pf.open( n_intref )
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()

    # compare optional output
    n_optout = 'rampfit1_opt_out_fitopt.fits'
    h = pf.open( n_optout )
    n_optref = _bigdata+'/miri/test_ramp_fit/jw00001001001_01101_00001_MIRIMAGE_opt.fits'
    href = pf.open( n_optref )
    newh = pf.HDUList([h['primary'],h['slope'],h['sigslope'],h['yint'],h['sigyint'],h['pedestal'],h['weights'],h['crmag']])
    newhref = pf.HDUList([href['primary'],href['slope'],href['sigslope'],href['yint'],href['sigyint'],href['pedestal'],href['weights'],href['crmag']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    assert result.identical, result.report()
