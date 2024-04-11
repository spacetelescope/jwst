"""Test ResidualFringeStep on MIRI MRS"""
import pytest

from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step


@pytest.mark.slow
@pytest.mark.bigdata
def test_residual_fringe_cal(rtdata, fitsdiff_default_kwargs):
    """Run residual fringe correction on MIRI IFUShort """

    input_file = 'jw01523001001_03101_00001_mirifushort_cal.fits'
    rtdata.get_data(f"miri/mrs/{input_file}")

    args = [
        'jwst.residual_fringe.ResidualFringeStep',
        input_file,
        '--save_results=true',
        '--skip=False'
    ]
    Step.from_cmdline(args)

    output = input_file.replace('cal', 'residual_fringe')
    rtdata.output = output

    rtdata.get_truth(f"truth/test_miri_residual_fringe/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
