import os

import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

OUTSTEM_BKG = "result_bkgasn"
OUTSTEM_SELFCAL = "result_selfcalasn"


@pytest.fixture(scope="module")
def run_pipeline_background(rtdata_module):
    """Run the pipeline with asn in which the background exposures are marked as `background`."""

    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01204-o021_20240127t024203_spec2_00010_asn.json")
    Step.from_cmdline(
        [
            "calwebb_spec2",
            rtdata.input,
            f"--steps.badpix_selfcal.output_file={OUTSTEM_BKG}",
            "--steps.badpix_selfcal.save_results=True",
            "--steps.badpix_selfcal.save_flagged_bkg=True",
            "--steps.badpix_selfcal.flagfrac_lower=0.0005",
            "--steps.badpix_selfcal.skip=False",
        ]
    )
    rtdata.output = f"{OUTSTEM_BKG}_badpix_selfcal.fits"
    return rtdata


@pytest.fixture(scope="module")
def run_pipeline_selfcal(rtdata_module):
    """Identical pipeline run to above, but input asn sets all background exposures as `selfcal` type."""
    rtdata = rtdata_module
    rtdata.get_asn("miri/mrs/jw01204-o021_20240127t024203_spec2_00010_selfcal_asn.json")
    Step.from_cmdline(
        [
            "calwebb_spec2",
            rtdata.input,
            f"--steps.badpix_selfcal.output_file={OUTSTEM_SELFCAL}",
            "--steps.badpix_selfcal.save_results=True",
            "--steps.badpix_selfcal.flagfrac_lower=0.0005",
            "--steps.badpix_selfcal.skip=False",
        ]
    )
    rtdata.output = f"{OUTSTEM_SELFCAL}_badpix_selfcal.fits"

    return rtdata


@pytest.mark.bigdata
def test_miri_mrs_badpix_selfcal(run_pipeline_selfcal, fitsdiff_default_kwargs):
    """Run a test for MIRI MRS data with dedicated background exposures."""

    rtdata = run_pipeline_selfcal

    # Get the truth file
    rtdata.get_truth(f"truth/test_miri_mrs_badpix_selfcal/{OUTSTEM_SELFCAL}_badpix_selfcal.fits")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    # check the bkg files in the background case, but not in the selfcal case
    for idx in range(4):
        fname = f"{OUTSTEM_SELFCAL}_badpix_selfcal_bkg_{idx}.fits"
        assert not os.path.isfile(fname)


@pytest.mark.parametrize(
    "basename",
    (
        [
            f"{OUTSTEM_BKG}_badpix_selfcal.fits",
        ]
        + [f"{OUTSTEM_BKG}_badpix_selfcal_bkg_{idx}.fits" for idx in range(4)]
    ),
)
@pytest.mark.bigdata
def test_miri_mrs_badpix_selfcal_bkg(basename, run_pipeline_background, fitsdiff_default_kwargs):
    """Run a test for MIRI MRS data with dedicated background exposures."""

    rtdata = run_pipeline_background

    # Get the truth file
    rtdata.output = basename
    rtdata.get_truth(f"truth/test_miri_mrs_badpix_selfcal/{basename}")

    # Compare the results and check the bkg files in the background case, but not in the selfcal case
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
