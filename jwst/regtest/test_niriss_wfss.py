"""Regression tests for NIRISS WFSS mode"""

import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_nis_wfss_spec2(rtdata_module, resource_tracker):
    """Run the calwebb_spec2 pipeline on NIRISS WFSS exposures"""
    rtdata = rtdata_module

    # These are the WFSS exposures we'll be processing
    spec2_asns = [
        "niriss/wfss/jw01324-o001_20220629t171902_spec2_001_asn.json",
        "niriss/wfss/jw01324-o001_20220629t171902_spec2_005_asn.json",
        "niriss/wfss/jw01324-o001_20220629t171902_spec2_007_asn.json",
    ]

    # Only the first exposure will be used to do comparisons
    rtdata.get_asn(spec2_asns[0])
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.bkg_subtract.save_results=true",
        "--steps.flat_field.save_results=true",
        "--steps.extract_2d.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.pathloss.save_results=true",
        "--steps.photom.save_results=true",
        "--steps.extract_1d.save_results=true",
        "--save_wfss_esec=true",
        "--steps.extract_2d.wfss_nbright=10",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    # Run the remaining exposures without doing comparisons, just so that
    # fresh results are available for level-3 processing
    for asn in spec2_asns[1:]:
        rtdata.get_asn(asn)
        args = ["calwebb_spec2", rtdata.input, "--steps.extract_2d.wfss_nbright=10"]
        Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_nis_wfss_spec3(run_nis_wfss_spec2, rtdata_module, resource_tracker):
    """Run the calwebb_spec3 pipeline"""
    rtdata = rtdata_module

    # Get the level3 association file and run the spec3 pipeline on it.
    # We don't need to retrieve any of the cal members of the association,
    # because they were all just created by the preceding spec2 test.
    rtdata.get_data("niriss/wfss/jw01324-o001_spec3_00005_asn.json")
    args = ["calwebb_spec3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_spec2(log_tracked_resources, run_nis_wfss_spec2):
    log_tracked_resources()


def test_log_tracked_resources_spec3(log_tracked_resources, run_nis_wfss_spec3):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix",
    ["assign_wcs", "bsub", "cal", "esec", "extract_2d", "flat_field", "photom", "srctype", "x1d"],
)
def test_nis_wfss_spec2(run_nis_wfss_spec2, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test for calwebb_spec2 applied to NIRISS WFSS data"""
    rtdata = rtdata_module
    rtdata.input = "jw01324001001_03101_00001_nis_rate.fits"
    output = "jw01324001001_03101_00001_nis_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_wfss/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("suffix", ["x1d", "c1d"])
def test_nis_wfss_spec3(run_nis_wfss_spec3, rtdata_module, suffix, fitsdiff_default_kwargs):
    """Regression test of the calwebb_spec3 pipeline applied to NIRISS WFSS data"""
    rtdata = rtdata_module
    rtdata.input = "jw01324-o001_spec3_00005_asn.json"
    output = "jw01324-o001_t001_niriss_f115w-gr150c_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_wfss/{output}")

    # Compare the results
    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
