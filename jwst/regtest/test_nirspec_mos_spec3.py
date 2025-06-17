import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
import numpy as np
from gwcs import wcstools

from jwst.stpipe import Step
from stdatamodels.jwst import datamodels

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module, resource_tracker):
    """Run calwebb_spec3 on NIRSpec MOS data."""
    rtdata = rtdata_module
    rtdata.get_asn("nirspec/mos/jw01345-o066_20230831t181155_spec3_00002_asn.json")

    # Run the calwebb_spec3 pipeline on the association
    args = ["calwebb_spec3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources(log_tracked_resources, run_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
@pytest.mark.parametrize(
    "source_id",
    [
        "b000000030",
        "b000000031",
        "s000004385",
        "s000007380",
        "v000000048",
        "v000000049",
        "v000000053",
        "v000000056",
    ],
)
def test_nirspec_mos_spec3(run_pipeline, suffix, source_id, fitsdiff_default_kwargs):
    """Check results of calwebb_spec3"""
    rtdata = run_pipeline

    output = f"jw01345-o066_{source_id}_nirspec_f170lp-g235m_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_mos_spec3/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "s2d":
        fitsdiff_default_kwargs["rtol"] = 1e-4
        fitsdiff_default_kwargs["atol"] = 1e-5

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    # Check output wavelength array against its own wcs
    if suffix == "s2d":
        tolerance = 1e-03
        dmr = datamodels.open(rtdata.output)

        w = dmr.meta.wcs
        x, y = wcstools.grid_from_bounding_box(w.bounding_box, step=(1, 1), center=True)
        _, _, wave = w(x, y)
        wlr = dmr.wavelength
        assert np.all(np.isclose(wave, wlr, atol=tolerance))
