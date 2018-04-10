import os
from astropy.io import fits as pf
from jwst.ipc.ipc_step import IPCStep

from ..helpers import add_suffix

BIGDATA = os.environ['TEST_BIGDATA']

def test_ipc_nircam():
    """Regression test of IPC step performed on NIRCam data."""

    output_file_base, output_file = add_suffix('ipc1_output.fits', 'ipc')

    try:
        os.remove(output_file)
    except OSError:
        pass

    IPCStep.call(BIGDATA + '/nircam/test_ipc_step/jw00017001001_01101_00001_NRCA3_uncal.fits',
                 output_file=output_file_base)
    h = pf.open(output_file)
    href = pf.open(BIGDATA + '/nircam/test_ipc_step/jw00017001001_01101_00001_NRCA3_ipc.fits')
    newh = pf.HDUList([h['primary'], h['sci']])
    newhref = pf.HDUList([href['primary'], href['sci']])
    result = pf.diff.FITSDiff(newh,
                                newhref,
                                ignore_keywords=['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                                rtol=1.e-6)
    h.close()
    href.close()
    result.report()
    try:
        assert result.identical
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)
