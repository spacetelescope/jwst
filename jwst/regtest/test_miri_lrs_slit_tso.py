import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]

SPEC2_BASENAME = "jw06219001001_04103_00001-seg001_mirimage"
ASN_ID = "o001"
TSO3_BASENAME = "jw06219-o001_t001_miri_p750l"


@pytest.fixture(scope="module")
def run_spec2_pipeline(rtdata_module, resource_tracker):
    """Run the calwebb_spec2 pipeline on a MIRI LRS FS TSO exposure."""
    rtdata = rtdata_module
    rtdata.get_data(f"miri/lrs/{SPEC2_BASENAME}_rateints.fits")

    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.flat_field.save_results=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso3_pipeline(run_spec2_pipeline, rtdata_module, resource_tracker):
    """Run the calwebb_tso3 pipeline on the output of run_spec2_pipeline."""
    rtdata = rtdata_module

    # Run spec2 for an extra input file
    second_file = SPEC2_BASENAME.replace("seg001", "seg002")
    rtdata.get_data(f"miri/lrs/{second_file}_rateints.fits")
    args = ["calwebb_spec2", rtdata.input]
    Step.from_cmdline(args)

    # Get the ASN for tso3
    rtdata.get_data("miri/lrs/jw06219-o001_20251211t073827_tso3_00001_asn_small.json")

    # Run tso3
    args = ["calwebb_tso3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_spec2(log_tracked_resources, run_spec2_pipeline):
    log_tracked_resources()


def test_log_tracked_resources_tso3(log_tracked_resources, run_tso3_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["assign_wcs", "srctype", "flat_field", "calints", "x1dints"])
def test_miri_lrs_slit_tso_spec2(
    run_spec2_pipeline, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Compare the output of a MIRI LRS FS TSO spec2 pipeline run."""
    rtdata = rtdata_module

    output_filename = f"{SPEC2_BASENAME}_{suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slit_tso/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_miri_lrs_slit_tso3_crfints(run_tso3_pipeline, rtdata_module, fitsdiff_default_kwargs):
    """Compare one of the crfints outputs from tso3 for MIRI LRS FS TSO."""
    rtdata = rtdata_module
    output_filename = f"{SPEC2_BASENAME}_{ASN_ID}_crfints.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slit_tso/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_miri_lrs_slit_tso3_x1dints(run_tso3_pipeline, rtdata_module, fitsdiff_default_kwargs):
    """Compare the x1dints output from tso3 for MIRI LRS FS TSO."""
    rtdata = rtdata_module

    output_filename = f"{TSO3_BASENAME}_x1dints.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slit_tso/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_miri_lrs_slit_tso3_whtlt(run_tso3_pipeline, rtdata_module, diff_astropy_tables):
    """Compare the whtlt output from tso3 for MIRI LRS FS TSO."""
    rtdata = rtdata_module

    output_filename = f"{TSO3_BASENAME}_whtlt.ecsv"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_miri_lrs_slit_tso/{output_filename}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)
