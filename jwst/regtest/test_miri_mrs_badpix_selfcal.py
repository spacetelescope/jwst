import pytest
from astropy.io.fits.diff import FITSDiff
from jwst.stpipe import Step
import os


@pytest.fixture(scope="module")
def run_pipeline_background(rtdata_module):

    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01204-o021_20240127t024203_spec2_00010_asn.json")
    Step.from_cmdline(['calwebb_spec2', rtdata.input,
                       "--steps.badpix_selfcal.save_results=True",
                       "--steps.badpix_selfcal.save_flagged_bkg=True",
                       "--steps.badpix_selfcal.flagfrac_lower=0.0005",
                       "--steps.badpix_selfcal.skip=False"])
    rtdata.output = "jw01204021001_02101_00004_mirifulong_badpix_selfcal.fits"
    return rtdata


@pytest.fixture(scope="module")
def run_pipeline_selfcal(rtdata_module):
    '''Identical pipeline run to above, but input asn sets all background exposures as `selfcal` type.
    '''
    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01204-o021_20240127t024203_spec2_00010_selfcal_asn.json")
    Step.from_cmdline(['calwebb_spec2', rtdata.input,
                       "--steps.badpix_selfcal.save_results=True",
                       "--steps.badpix_selfcal.flagfrac_lower=0.0005",
                       "--steps.badpix_selfcal.skip=False"])
    rtdata.output = "jw01204021001_02101_00004_mirifulong_badpix_selfcal.fits"

    return rtdata


@pytest.mark.bigdata
def test_miri_mrs_badpix_selfcal(run_pipeline_selfcal, fitsdiff_default_kwargs):
    """Run a test for MIRI MRS data with dedicated background exposures."""

    rtdata = run_pipeline_selfcal

    # Get the truth file
    rtdata.get_truth("truth/test_miri_mrs_badpix_selfcal/jw01204021001_02101_00004_mirifulong_badpix_selfcal.fits")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    # check the bkg files in the background case, but not in the selfcal case
    for idx in range(4):
        fname = f"jw01204021001_02101_00004_mirifulong_badpix_selfcal_bkg_{idx}.fits"
        assert not os.path.isfile(fname)


@pytest.mark.bigdata
def test_miri_mrs_badpix_selfcal_bkg(run_pipeline_background, fitsdiff_default_kwargs):
    """Run a test for MIRI MRS data with dedicated background exposures."""

    rtdata = run_pipeline_background

    # Get the truth file
    rtdata.get_truth("truth/test_miri_mrs_badpix_selfcal/jw01204021001_02101_00004_mirifulong_badpix_selfcal.fits")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    # check the bkg files in the background case, but not in the selfcal case
    for idx in range(4):
        fname = f"jw01204021001_02101_00004_mirifulong_badpix_selfcal_bkg_{idx}.fits"
        truth = rtdata.get_truth(f"truth/test_miri_mrs_badpix_selfcal/{fname}")
        diff = FITSDiff(fname, truth, **fitsdiff_default_kwargs)
        assert diff.identical, diff.report()