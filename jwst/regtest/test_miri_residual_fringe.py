"""Test ResidualFringeStep on MIRI MRS"""
import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope='module')
def run_residual_fringe_output(rtdata_module):
    """Run residual_fringe t"""
    rtdata = rtdata_module
    rtdata.get_asn('miri/mrs/rfc_asn.json')

    args = [
        'jwst.residual_fringe.ResidualFringeStep',
        rtdata.input,
        '--save_results=true'
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.slow
@pytest.mark.parametrize(
    'output',
    ['outtest_residual_fringe.fits']
)
def test(run_cube_build_single_output, output, fitsdiff_default_kwargs):
    """Test just running residual fringe correction"""
    rtdata = run_residual_fringe_output
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(f'truth/test_residual_fringe/{output}')

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


#@pytest.mark.bigdata
#def test_residual_fringe_cal(rtdata, fitsdiff_default_kwargs):
#    """Run residual fringe correction on MIRI IFUShort
#    input_file = 'det_image_seq2_MIRIFUSHORT_12SHORTexp1_cal.fits'
#    rtdata.get_data(f"miri/mrs/{input_file}")

#    args = [
#        'jwst.residual_fringe.ResidualFringeStep',
#        input_file,
#        '--save_results=true',
#    ]
#    Step.from_cmdline(args)

#    output = input_file.replace('cal', 'residual_fringe')
#    rtdata.output = output

#    rtdata.get_truth(f"truth/test_miri_residual_fringe/{output}")

#    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
#    assert diff.identical, diff.report()
