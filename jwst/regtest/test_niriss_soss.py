import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_tso_spec2(jail, rtdata_module):
    """Run stage 2 pipeline on NIRISS SOSS data."""
    rtdata = rtdata_module
    collect_pipeline_cfgs("config")

    # Run tso-spec2 pipeline on the first _rateints file, saving intermediate products
    rtdata.get_data("niriss/soss/jw00625023001_03101_00001-seg001_nis_rateints.fits")
    args = ["config/calwebb_tso-spec2.cfg", rtdata.input,
        "--steps.flat_field.save_results=True",
        "--steps.srctype.save_results=True",
        ]
    Step.from_cmdline(args)

    # Run tso-spec2 pipeline on the second _rateints file, without saving or
    # checking any results (simply create a fresh input for level-3 test)
    rtdata.get_data("niriss/soss/jw00625023001_03101_00001-seg002_nis_rateints.fits")
    args = ["config/calwebb_tso-spec2.cfg", rtdata.input]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_tso_spec3(jail, rtdata_module, run_tso_spec2):
    """Run stage 3 pipeline on NIRISS SOSS data."""
    rtdata = rtdata_module
    # Get the level3 assocation json file (though not its members) and run
    # the tso3 pipeline on all _calints files listed in association
    rtdata.get_data("niriss/soss/jw00625-o023_20191210t204036_tso3_001_asn.json")
    args = ["config/calwebb_tso3.cfg", rtdata.input]
    Step.from_cmdline(args)

@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["calints", "flat_field", "srctype", "x1dints"])
def test_niriss_soss_stage2(rtdata_module, run_tso_spec2, fitsdiff_default_kwargs, suffix):
    """Regression test of tso-spec2 pipeline performed on NIRISS SOSS data."""
    rtdata = rtdata_module
    output = f"jw00625023001_03101_00001-seg001_nis_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_niriss_soss_stage3_crfints(rtdata_module, run_tso_spec3, fitsdiff_default_kwargs):
    """Regression test of tso-spec3 pipeline outlier_detection results performed on NIRISS SOSS data."""
    rtdata = rtdata_module
    output = "jw00625023001_03101_00001-seg001_nis_o023_crfints.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_niriss_soss_stage3_x1dints(run_tso_spec3, rtdata_module, fitsdiff_default_kwargs):
    """Regression test of tso-spec3 pipeline extract_1d results performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = "jw00625-o023_t001_niriss_clear-gr700xd-substrip256_x1dints.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_niriss_soss_stage3_whtlt(run_tso_spec3, rtdata_module, diff_astropy_tables):
    """Regression test of tso-spec3 pipeline white_light results performed on NIRISS SOSS data."""
    rtdata = rtdata_module

    output = "jw00625-o023_t001_niriss_clear-gr700xd-substrip256_whtlt.ecsv"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_soss_stages/{output}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth)
