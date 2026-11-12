"""Regression tests for NIRISS WFSS with WCS tweaks."""

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_nis_wfss_direct_image3(rtdata_module, resource_tracker):
    """Run the calwebb_image3 pipeline on NIRISS WFSS direct images."""
    rtdata = rtdata_module

    image3_asn = "niriss/wfss/jw01510-o001_20261017t073409_image3_00001_asn.json"
    rtdata.get_asn(image3_asn)

    # Run image3 on a set of images including one visit with a failed
    # guide star acquisition: visit 4 has a significant tweak applied.
    # Set bright=120 to make sure alignment succeeds.
    args = [
        "calwebb_image3",
        rtdata.input,
        "--steps.tweakreg.bright=120",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_nis_wfss_spec2(run_nis_wfss_direct_image3, rtdata_module, resource_tracker):
    """Run the calwebb_spec2 pipeline on NIRISS WFSS exposures needing a WCS tweak."""
    rtdata = rtdata_module

    # Get the association + rate file only, use image3 input as created here
    spec2_rate = "niriss/wfss/jw01510001004_0210h_00001_nis_rate.fits"
    spec2_asn = "niriss/wfss/jw01510-o001_20261017t073409_spec2_00001_asn.json"
    rtdata.get_data(spec2_rate)
    rtdata.get_data(spec2_asn)

    # Run spec3 extracting one particular source: it should be clearly centered if
    # tweak is correctly applied
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.extract_2d.save_results=true",
        "--steps.extract_2d.source_ra=81.08343554630416",
        "--steps.extract_2d.source_dec=-70.10274149882585",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.mark.parametrize("suffix", ["i2d", "segm"])
def test_nis_wfss_direct_tweakreg(
    run_nis_wfss_direct_image3, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Check output products for NIRISS WFSS direct images."""
    rtdata = rtdata_module
    output = f"jw01510-o001_t001_niriss_clear-f200w_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_wfss_tweakreg/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nis_wfss_direct_tweakreg_catalog(
    run_nis_wfss_direct_image3, rtdata_module, diff_astropy_tables
):
    """Check output catalog for NIRISS WFSS direct images."""
    rtdata = rtdata_module
    output = "jw01510-o001_t001_niriss_clear-f200w_cat.ecsv"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_wfss_tweakreg/{output}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth, rtol=1e-4, atol=1e-5)


@pytest.mark.parametrize("suffix", ["extract_2d", "cal", "x1d"])
def test_nis_wfss_spec2_tweakreg(
    run_nis_wfss_spec2, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Check output products for NIRISS WFSS data with a significant tweak to WCS."""
    rtdata = rtdata_module
    output = f"jw01510001004_0210h_00001_nis_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_wfss_tweakreg/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
