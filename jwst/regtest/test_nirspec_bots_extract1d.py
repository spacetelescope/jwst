import os
import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_extract(rtdata_module, request):
    """Run the extract1d step with custom parameters on BOTS data (S1600A1 slit)."""

    rtdata = rtdata_module

    # Get the custom reference file and input exposure
    ref_file = 'jwst_nirspec_extract1d_custom_g395h.json'
    rtdata.get_data(f'nirspec/tso/{ref_file}')
    rtdata.get_data('nirspec/tso/jw01118005001_04101_00001-first20_nrs1_calints.fits')

    # Run the calwebb_spec2 pipeline;
    args = ["extract_1d", rtdata.input,
            f"--override_extract1d={ref_file}",
            "--suffix=x1dints"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
def test_nirspec_bots_custom_extraction(run_extract, fitsdiff_default_kwargs):
    """
        Regression test of calwebb_spec2 pipeline performed on NIRSpec
        fixed-slit data that uses the NRS_BRIGHTOBJ mode (S1600A1 slit).
    """
    rtdata = run_extract
    output = 'jw01118005001_04101_00001-first20_nrs1_x1dints.fits'
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_bots_extract1d", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
