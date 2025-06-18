import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module, resource_tracker):
    """Run calwebb_spec3 on NIRSpec MOS data."""
    rtdata = rtdata_module
    rtdata.get_asn("nirspec/mos/jw02674-o004_20240305t054741_spec3_00001_asn.json")

    # Run the calwebb_spec3 pipeline on the association
    args = ["calwebb_spec3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_spec3(log_tracked_resources, run_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
@pytest.mark.parametrize(
    "source_id",
    [
        "b000000003",
        "b000000004",
        "b000000048",
        "b000000052",
        "s000001354",
        "s000012105",
        "s000034946",
    ],
)
def test_nirspec_mos_fs_spec3(run_pipeline, suffix, source_id, fitsdiff_default_kwargs):
    """Check results of calwebb_spec3"""
    rtdata = run_pipeline

    output = f"jw02674-o004_{source_id}_nirspec_f290lp-g395m_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_mos_fs_spec3/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "s2d":
        fitsdiff_default_kwargs["rtol"] = 1e-4
        fitsdiff_default_kwargs["atol"] = 1e-5

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
