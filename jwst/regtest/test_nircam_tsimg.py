import pytest
from astropy.io.fits.diff import FITSDiff
from astropy.table import Table, setdiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipelines(jail, rtdata_module):
    """Run stage 2 and 3 pipelines on NIRCam TSIMG data."""

    rtdata = rtdata_module

    collect_pipeline_cfgs("config")

    # Run the calwebb_tso-image2 pipeline on each of the 2 inputs
    rate_files = [
        "nircam/tsimg/jw00312006001_02102_00001-seg001_nrcb1_rateints.fits",
        "nircam/tsimg/jw00312006001_02102_00001-seg002_nrcb1_rateints.fits"
    ]
    for rate_file in rate_files:
        rtdata.get_data(rate_file)
        args = ["config/calwebb_tso-image2.cfg", rtdata.input]
        Step.from_cmdline(args)

    # Get the level3 assocation json file (though not its members) and run
    # the tso3 pipeline on all _calints files listed in association
    rtdata.get_data("nircam/tsimg/jw00312-o006_20191225t115310_tso3_001_asn.json")
    args = ["config/calwebb_tso3.cfg", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["calints", "o006_crfints"])
def test_nircam_tsimg_stage2(run_pipelines, fitsdiff_default_kwargs, suffix):
    """Regression test of tso-image2 pipeline performed on NIRCam TSIMG data."""
    rtdata = run_pipelines
    rtdata.input = "jw00312006001_02102_00001-seg001_nrcb1_rateints.fits"
    output = "jw00312006001_02102_00001-seg001_nrcb1_" + suffix + ".fits"
    rtdata.output = output

    rtdata.get_truth("truth/test_nircam_tsimg_stage23/" + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_tsimage_stage3_phot(run_pipelines):
    rtdata = run_pipelines
    rtdata.input = "jw00312-o006_20191225t115310_tso3_001_asn.json"
    rtdata.output = "jw00312-o006_t001_nircam_f210m-clear-sub64p_phot.ecsv"
    rtdata.get_truth("truth/test_nircam_tsimg_stage23/jw00312-o006_t001_nircam_f210m-clear-sub64p_phot.ecsv")

    table = Table.read(rtdata.output)
    table_truth = Table.read(rtdata.truth)

    # setdiff returns a table of length zero if there is no difference
    assert len(setdiff(table, table_truth)) == 0
