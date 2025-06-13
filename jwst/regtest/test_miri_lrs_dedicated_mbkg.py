import os
import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.master_background import MasterBackgroundStep


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    rtdata = rtdata_module

    rtdata.get_asn("miri/lrs/jw01529-o003_spec3_with_bg_00001_asn.json")

    MasterBackgroundStep.call(rtdata.input, save_results=True, suffix="mbsub", save_background=True)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "output",
    [
        "jw01529-o004_t002_miri_p750l_o003_masterbg1d.fits",
        "jw01529003001_03103_00011_mirimage_o003_masterbg2d.fits",
        "jw01529003001_03103_00011_mirimage_mbsub.fits",
    ],
)
def test_miri_lrs_dedicated_mbkg(run_pipeline, fitsdiff_default_kwargs, output):
    """Run a test for MIRI LRS data with dedicated background exposures."""

    rtdata = run_pipeline
    rtdata.output = output

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_miri_lrs_dedicated_mbkg", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
