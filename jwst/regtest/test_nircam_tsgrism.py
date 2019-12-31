import pytest
from astropy.io.fits.diff import FITSDiff
from astropy.table import Table, setdiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipelines(jail, rtdata_module):
    """Run stage 2-3 tso pipelines on NIRCAM TSO grism data."""
    rtdata = rtdata_module
    collect_pipeline_cfgs("config")

    # Run tso-spec2 pipeline on the _rateints file, saving intermediate products
    rtdata.get_data("nircam/tsgrism/jw00721012001_03103_00001-seg001_nrcalong_rateints.fits")
    args = ["config/calwebb_tso-spec2.cfg", rtdata.input,
        "--steps.flat_field.save_results=True",
        "--steps.extract_2d.save_results=True",
        "--steps.srctype.save_results=True",
        ]
    Step.from_cmdline(args)

    # Get the level3 assocation json file (though not its members) and run
    # the tso3 pipeline on all _calints files listed in association
    rtdata.get_data("nircam/tsgrism/jw00721-o012_20191119t043909_tso3_001_asn.json")
    args = ["config/calwebb_tso3.cfg", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["calints", "extract_2d", "flat_field",
    "o012_crfints", "srctype", "x1dints"])
def test_nircam_tsgrism_stage2(run_pipelines, fitsdiff_default_kwargs, suffix):
    """Regression test of tso-spec2 pipeline performed on NIRCam TSO grism data."""
    rtdata = run_pipelines
    rtdata.input = "jw00721012001_03103_00001-seg001_nrcalong_rateints.fits"
    output = "jw00721012001_03103_00001-seg001_nrcalong_" + suffix + ".fits"
    rtdata.output = output

    rtdata.get_truth("truth/test_nircam_tsgrism_stages/" + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_tsgrism_stage3_x1dints(run_pipelines, fitsdiff_default_kwargs):
    rtdata = run_pipelines
    rtdata.input = "jw00721-o012_20191119t043909_tso3_001_asn.json"
    rtdata.output = "jw00721-o012_t004_nircam_f444w-grismr-subgrism256_x1dints.fits"
    rtdata.get_truth("truth/test_nircam_tsgrism_stages/jw00721-o012_t004_nircam_f444w-grismr-subgrism256_x1dints.fits")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_tsgrism_stage3_whtlt(run_pipelines):
    rtdata = run_pipelines
    rtdata.input = "jw00721-o012_20191119t043909_tso3_001_asn.json"
    rtdata.output = "jw00721-o012_t004_nircam_f444w-grismr-subgrism256_whtlt.ecsv"
    rtdata.get_truth("truth/test_nircam_tsgrism_stages/jw00721-o012_t004_nircam_f444w-grismr-subgrism256_whtlt.ecsv")

    table = Table.read(rtdata.output)
    table_truth = Table.read(rtdata.truth)

    # setdiff returns a table of length zero if there is no difference
    assert len(setdiff(table, table_truth)) == 0
