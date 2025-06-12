import os
import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.master_background import MasterBackgroundStep


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    rtdata = rtdata_module

    rtdata.get_asn("miri/mrs/jw01031-c1004_20241028t205539_spec3_subset_asn.json")

    MasterBackgroundStep.call(rtdata.input, save_results=True, suffix="mbsub", save_background=True)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("exposure", [1, 2, 3, 4])
@pytest.mark.parametrize("suffix", ["mbsub", "c1004_masterbg2d"])
def test_miri_mrs_dedicated_mbkg(run_pipeline, fitsdiff_default_kwargs, exposure, suffix):
    """Run a test for MIRI MRS data with dedicated background exposures."""
    rtdata = run_pipeline
    output_file = f"jw01031001001_02101_0000{exposure}_mirifulong_{suffix}.fits"
    rtdata.output = output_file

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_miri_mrs_dedicated_mbkg", output_file))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_mrs_dedicated_masterbg1d(run_pipeline, fitsdiff_default_kwargs):
    rtdata = run_pipeline

    # Check 1D masterbg output product: created with root name from nod 1
    output_file = "jw01031006001_02101_00001_mirifulong_c1004_masterbg1d.fits"
    rtdata.output = output_file

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_miri_mrs_dedicated_mbkg", output_file))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
