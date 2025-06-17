import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module, resource_tracker):
    """Run the calwebb_spec2 pipeline on a single NIRSpec MOS/FS exposure."""

    rtdata = rtdata_module

    # Get the MSA metadata file referenced in the input exposure
    rtdata.get_data("nirspec/mos/jw02674004001_01_msa.fits")

    # Get the input ASN file and exposures
    rtdata.get_data("nirspec/mos/jw02674004001_03101_00001_nrs1_rate.fits")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.msa_flagging.save_results=true",
        "--steps.master_background_mos.save_results=true",
        "--steps.extract_2d.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.wavecorr.save_results=true",
        "--steps.flat_field.save_results=true",
        "--steps.flat_field.save_interpolated_flat=true",
        "--steps.pathloss.save_results=true",
        "--steps.barshadow.save_results=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_spec2(log_tracked_resources, run_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix",
    [
        "assign_wcs",
        "msa_flagging",
        "extract_2d",
        "srctype",
        "master_background_mos",
        "wavecorr",
        "flat_field",
        "interpolatedflat",
        "pathloss",
        "barshadow",
        "wavecorr_fs",
        "flat_field_fs",
        "interpolatedflat_fs",
        "pathloss_fs",
        "cal",
        "s2d",
        "x1d",
    ],
)
def test_nirspec_mos_fs_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test for calwebb_spec2 on a NIRSpec MOS/FS exposure."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = f"jw02674004001_03101_00001_nrs1_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nirspec_mos_fs_spec2/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
