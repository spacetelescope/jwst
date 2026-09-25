import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_pipelines(rtdata_module, resource_tracker):
    """Run stage 2 and 3 pipelines on MIRI TSO image data."""

    rtdata = rtdata_module

    # Run the calwebb_image2 pipeline on each of the 2 inputs
    rate_files = [
        "miri/image/jw01177007001_03101_00001-seg001_mirimage_rateints.fits",
        "miri/image/jw01177007001_03101_00001-seg002_mirimage_rateints.fits",
    ]
    for rate_file in rate_files:
        rtdata.get_data(rate_file)
        # With default values, assign_wcs raises a warning for these data,
        # related to SIP approximation accuracy. Decrease the accuracy
        # for these tests.
        sip_max_error = 0.05
        args = [
            "calwebb_image2",
            rtdata.input,
            f"--steps.assign_wcs.sip_max_pix_error={sip_max_error}",
            f"--steps.assign_wcs.sip_max_inv_pix_error={sip_max_error}",
        ]

        Step.from_cmdline(args)

    # Get the level3 association json file (though not its members) and run
    # the tso3 pipeline on all _calints files listed in association
    rtdata.get_data("miri/image/jw01177-o007_tso3_00001_asn.json")
    args = ["calwebb_tso3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_miri_image_tso(log_tracked_resources, run_pipelines):
    log_tracked_resources()


@pytest.mark.parametrize("segment", ["seg001", "seg002"])
@pytest.mark.parametrize("suffix", ["calints", "o007_crfints"])
def test_miri_image_tso_exposure_data(run_pipelines, fitsdiff_default_kwargs, segment, suffix):
    """Regression test of tso-image2 pipeline performed on MIRI image TSO data."""
    rtdata = run_pipelines
    rtdata.input = f"jw01177007001_03101_00001-{segment}_mirimage_rateints.fits"
    output = f"jw01177007001_03101_00001-{segment}_mirimage_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth("truth/test_miri_image_tso/" + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_miri_image_tso_stage3_phot(run_pipelines, diff_astropy_tables):
    rtdata = run_pipelines
    rtdata.input = "jw01177-o007_tso3_00001_asn.json"
    rtdata.output = "jw01177-o007_t005_miri_f1500w_phot.ecsv"
    rtdata.get_truth("truth/test_miri_image_tso/jw01177-o007_t005_miri_f1500w_phot.ecsv")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)
