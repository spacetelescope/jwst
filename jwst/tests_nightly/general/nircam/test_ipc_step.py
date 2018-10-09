import pytest
from astropy.io import fits
from jwst.ipc.ipc_step import IPCStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_ipc_nircam(_bigdata):
    """
    Regression test of IPC step performed on NIRCam data.
    """
    suffix = 'ipc'
    output_file_base, output_file = add_suffix('ipc1_output.fits', suffix)

    IPCStep.call(_bigdata + '/nircam/test_ipc_step/jw00017001001_01101_00001_NRCA3_uncal.fits',
                 output_file=output_file_base, suffix=suffix
                 )
    h = fits.open(output_file)
    href = fits.open(_bigdata + '/nircam/test_ipc_step/jw00017001001_01101_00001_NRCA3_ipc.fits')
    newh = fits.HDUList([h['primary'], h['sci']])
    newhref = fits.HDUList([href['primary'], href['sci']])
    result = fits.diff.FITSDiff(newh,
                                newhref,
                                ignore_keywords=['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                                rtol=1.e-6)
    assert result.identical, result.report()
