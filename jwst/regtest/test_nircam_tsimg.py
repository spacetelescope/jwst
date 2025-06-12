import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.lib.set_telescope_pointing import add_wcs
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_pipelines(rtdata_module, resource_tracker):
    """Run stage 2 and 3 pipelines on NIRCam TSO image data."""

    rtdata = rtdata_module

    # Run the calwebb_tso-image2 pipeline on each of the 2 inputs
    rate_files = [
        "nircam/tsimg/jw01068006001_03103_00001-seg001_nrcb1_rateints.fits",
    ]
    for rate_file in rate_files:
        rtdata.get_data(rate_file)
        args = ["calwebb_image2", rtdata.input]
        Step.from_cmdline(args)

    # Get the level3 association json file (though not its members) and run
    # the tso3 pipeline on all _calints files listed in association
    rtdata.get_data("nircam/tsimg/jw01068-o006_20240401t151322_tso3_00002_asn.json")
    args = ["calwebb_tso3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_tsimg(log_tracked_resources, run_pipelines):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["calints", "o006_crfints"])
def test_nircam_tsimg_stage2(run_pipelines, fitsdiff_default_kwargs, suffix):
    """Regression test of tso-image2 pipeline performed on NIRCam TSIMG data."""
    rtdata = run_pipelines
    rtdata.input = "jw01068006001_03103_00001-seg001_nrcb1_rateints.fits"
    output = "jw01068006001_03103_00001-seg001_nrcb1_" + suffix + ".fits"
    rtdata.output = output

    rtdata.get_truth("truth/test_nircam_tsimg_stage23/" + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nircam_tsimage_stage3_phot(run_pipelines, diff_astropy_tables):
    rtdata = run_pipelines
    rtdata.input = "jw01068-o006_20240401t151322_tso3_00002_asn.json"
    rtdata.output = "jw01068-o006_t004_nircam_f150w2-f164n-sub64p_phot.ecsv"
    rtdata.get_truth(
        "truth/test_nircam_tsimg_stage23/jw01068-o006_t004_nircam_f150w2-f164n-sub64p_phot.ecsv"
    )

    assert diff_astropy_tables(rtdata.output, rtdata.truth)


def test_nircam_setpointing_tsimg(rtdata, fitsdiff_default_kwargs):
    """
    Regression test of the set_telescope_pointing script on a level-1b
    NIRCam TSO imaging file.
    """

    rtdata.get_data("nircam/tsimg/jw06780001001_02103_00001-seg002_nrcblong_uncal.fits")
    # The add_wcs function overwrites its input, so output = input
    rtdata.output = rtdata.input

    # Call the WCS routine
    add_wcs(rtdata.input)

    rtdata.get_truth(
        "truth/test_nircam_setpointing/jw06780001001_02103_00001-seg002_nrcblong_uncal.fits"
    )

    fitsdiff_default_kwargs["rtol"] = 1e-6
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
