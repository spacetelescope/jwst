import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_detector1pipeline(rtdata_module, resource_tracker):
    """Run calwebb_detector1 on NIRCam imaging long data"""
    rtdata = rtdata_module
    rtdata.get_data("nircam/dhs/jw04453010001_02106_00001_nrca1_genheader_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files
    args = [
        "calwebb_detector1",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.superbias.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.jump.save_results=True",
        "--steps.jump.rejection_threshold=50.0",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_det1(log_tracked_resources, run_detector1pipeline):
    log_tracked_resources()


@pytest.fixture(scope="module")
def run_spec2pipeline(rtdata_module, resource_tracker):
    """Run calwebb_detector1 on NIRCam imaging long data"""
    rtdata = rtdata_module
    rtdata.get_data("nircam/dhs/tpauly_nrca1_regions.asdf")
    rtdata.get_data("nircam/dhs/nircam_nrca1_specwcs.asdf")
    rtdata.get_data("nircam/dhs/jw04453010001_02106_00001_nrca1_genheader_rateints.fits")

    # Run spec2 pipeline on rateints file
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=True",
        "--steps.assign_wcs.override_regions=tpauly_nrca1_regions.asdf",
        "--steps.assign_wcs.override_specwcs=nircam_nrca1_specwcs.asdf",
        "--steps.extract_2d.save_results=True",
        "--steps.srctype.save_results=True",
        "--steps.photom.save_results=True",
        "--steps.extract_1d.save_results=True",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_spec2(log_tracked_resources, run_spec2pipeline):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix",
    [
        "dq_init",
        "saturation",
        "superbias",
        "refpix",
        "linearity",
        "dark_current",
        "jump",
        "rate",
        "rateints",
    ],
)
def test_nircam_dhs_detector1(
    run_detector1pipeline, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Regression test of detector1 pipeline performed on NIRCam DHS mock data."""
    rtdata = rtdata_module
    rtdata.input = "jw04453010001_02106_00001_nrca1_genheader_uncal.fits"
    output = "jw04453010001_02106_00001_nrca1_genheader_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_dhs/{output}")

    # Ignore the custom dark because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_DARK")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "suffix",
    [
        "assign_wcs",
        "extract_2d",
        "srctype",
        "photom",
        "calints",
        "x1dints",
    ],
)
def test_nircam_dhs_spec2(run_spec2pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 pipeline performed on NIRCam DHS mock data."""
    rtdata = rtdata_module
    rtdata.input = "jw04453010001_02106_00001_nrca1_genheader_rateints.fits"
    output = "jw04453010001_02106_00001_nrca1_genheader_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_dhs/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
