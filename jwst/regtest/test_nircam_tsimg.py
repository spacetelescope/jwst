import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipelines(rtdata_module):
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
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
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


@pytest.mark.bigdata
def test_nircam_tsimage_stage3_phot(run_pipelines, diff_astropy_tables):
    rtdata = run_pipelines
    rtdata.input = "jw01068-o006_20240401t151322_tso3_00002_asn.json"
    rtdata.output = "jw01068-o006_t004_nircam_f150w2-f164n-sub64p_phot.ecsv"
    rtdata.get_truth("truth/test_nircam_tsimg_stage23/jw01068-o006_t004_nircam_f150w2-f164n-sub64p_phot.ecsv")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)
