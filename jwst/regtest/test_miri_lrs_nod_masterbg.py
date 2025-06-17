import os
import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline_with_master_bg(rtdata_module):
    rtdata = rtdata_module

    # Get rate files and custom extraction parameters
    rate_files = [
        "jw01530005001_03103_00001_mirimage_rate.fits",
        "jw01530005001_03103_00002_mirimage_rate.fits",
    ]
    ref_files = ["jwst_miri_extract1d_nod1_bg.json", "jwst_miri_extract1d_nod2_bg.json"]
    for rate, ref in zip(rate_files, ref_files):
        rtdata.get_data(f"miri/lrs/{ref}")
        rtdata.get_data(f"miri/lrs/{rate}")

        Step.from_cmdline(
            [
                "calwebb_spec2",
                rate,
                f"--steps.extract_1d.override_extract1d={ref}",
                "--steps.extract_1d.use_source_pos=False",
            ]
        )

    # Get the ASN, but not the data - use the spec2 input products
    rtdata.get_data("miri/lrs/jw01530-o005_20221202t204827_spec3_00001_asn.json")
    Step.from_cmdline(
        [
            "calwebb_spec3",
            rtdata.input,
            "--steps.master_background.skip=False",
            "--steps.master_background.save_background=True",
            "--steps.master_background.save_results=True",
        ]
    )

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("nod", [1, 2])
@pytest.mark.parametrize("suffix", ["cal", "x1d", "s2d", "mbsub", "o005_crf", "o005_masterbg2d"])
def test_miri_lrs_nod_bg(run_pipeline_with_master_bg, fitsdiff_default_kwargs, nod, suffix):
    """Run a regression test for nodded MIRI LRS data with background extraction."""
    rtdata = run_pipeline_with_master_bg

    # Check output products for nod 1 or 2
    output_file = f"jw01530005001_03103_0000{nod}_mirimage_{suffix}.fits"
    rtdata.output = output_file

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_miri_lrs_nod_masterbg", output_file))

    # Compare the results
    # Ignore the custom extract1d file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_EXTR1D")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_lrs_nod_masterbg1d(run_pipeline_with_master_bg, fitsdiff_default_kwargs):
    rtdata = run_pipeline_with_master_bg

    # Check 1D masterbg output product: created with root name from nod 1
    output_file = "jw01530005001_03103_00001_mirimage_o005_masterbg1d.fits"
    rtdata.output = output_file

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_miri_lrs_nod_masterbg", output_file))

    # Compare the results
    # Ignore the custom extract1d file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_EXTR1D")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["s2d", "x1d"])
def test_miri_lrs_nod_bg_spec3(run_pipeline_with_master_bg, fitsdiff_default_kwargs, suffix):
    """Run a regression test for nodded MIRI LRS with master background subtraction."""
    rtdata = run_pipeline_with_master_bg

    # Check combined output products
    output_file = f"jw01530-o005_t004_miri_p750l_{suffix}.fits"
    rtdata.output = output_file

    # Get the truth file
    rtdata.get_truth(os.path.join("truth/test_miri_lrs_nod_masterbg", output_file))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
