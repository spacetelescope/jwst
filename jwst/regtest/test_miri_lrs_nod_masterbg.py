import os
import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.master_background import MasterBackgroundStep


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):

    rtdata = rtdata_module

    rtdata.get_asn("miri/lrs/miri_lrs_mbkg_nodded_spec3_asn.json")

    MasterBackgroundStep.call(rtdata.input, save_results=True, suffix='master_background')

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("nod_seq", ["seq1", "seq2"])
def test_miri_lrs_nod_masterbg(run_pipeline, fitsdiff_default_kwargs, nod_seq):
    """Run a regression test for nodded MIRI LRS data."""

    rtdata = run_pipeline
    rtdata.output = "miri_lrs_nod_" + nod_seq + "_exp1_master_background.fits"

    # Get the truth file
    rtdata.get_truth(os.path.join(
        "truth/test_miri_lrs_nod_masterbg",
        "miri_lrs_nod_" + nod_seq + "_exp1_master_background.fits"))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
