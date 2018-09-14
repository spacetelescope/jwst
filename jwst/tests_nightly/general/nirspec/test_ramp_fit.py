import pytest
from astropy.io import fits
from jwst.ramp_fitting.ramp_fit_step import RamfitsitStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_ramp_fit_nirspec(_bigdata):
    """

    Regression test of ramp_fit step performed on NIRSpec data. This is a single
    integration dataset.

    """
    output_file_base, output_files = add_suffix('ramfitsit_output.fits', 'ramfitsit', list(range(1)))

    RamfitsitStep.call(_bigdata+'/nirspec/test_ramp_fit/jw00023001001_01101_00001_NRS1_jump.fits',
                      output_file=output_file_base,
                      save_opt=True,
                      opt_name='ramfitsit_opt_out.fits', name='Ramfitsit'
                      )

    # compare primary output
    n_priout = output_files[0]
    h = fits.open( n_priout )
    n_priref = _bigdata+'/nirspec/test_ramp_fit/jw00023001001_01101_00001_NRS1_ramp_fit.fits'
    href = fits.open( n_priref )
    newh = fits.HDUList([h['primary'],h['sci'],h['err'],h['dq']])
    newhref = fits.HDUList([href['primary'],href['sci'],href['err'],href['dq']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )


    assert result.identical, result.report()


    # compare optional output
    n_optout = 'ramfitsit_opt_out_fitopt.fits'
    h = fits.open( n_optout )
    n_optref = _bigdata+'/nirspec/test_ramp_fit/jw00023001001_01101_00001_NRS1_opt.fits'
    href = fits.open( n_optref )
    newh = fits.HDUList([h['primary'],h['slope'],h['sigslope'],h['yint'],h['sigyint'],h['pedestal'],h['weights'],h['crmag']])
    newhref = fits.HDUList([href['primary'],href['slope'],href['sigslope'],href['yint'],href['sigyint'],href['pedestal'],href['weights'],href['crmag']])
    result = fits.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )


    assert result.identical, result.report()
