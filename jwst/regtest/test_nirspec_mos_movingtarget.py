import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_spec2_pipeline(rtdata_module):
    """Run the calwebb_spec2 pipeline on a single NIRSpec MOS Moving Target exposure."""

    rtdata = rtdata_module

    # Get the MSA metadata file referenced in the input exposure
    rtdata.get_data("nirspec/mos/jw01444003001_01_msa.fits")

    # Get the input ASN file and exposures
    rtdata.get_data("nirspec/mos/jw01444003001_04102_00001_nrs1_rate.fits")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.msa_flagging.save_results=true",
        "--steps.master_background_mos.save_results=true",
        "--steps.extract_2d.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.wavecorr.save_results=true",
        "--steps.flat_field.save_results=true",
        "--steps.pathloss.save_results=true",
        "--steps.barshadow.save_results=true",
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.fixture(scope="module")
def run_spec3_pipeline(run_spec2_pipeline, rtdata_module):
    """Run calwebb_spec3 on NIRSpec MOS data."""
    rtdata = rtdata_module
    rtdata.get_data("nirspec/mos/jw01444-o003_20240815t170722_spec3_00001_asn.json")

    # Run the calwebb_spec3 pipeline on the association
    args = ["calwebb_spec3", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "suffix",
    [
        "assign_wcs",
        "msa_flagging",
        "extract_2d",
        "srctype",
        "master_background_mos",
        "wavecorr",
        "flat_field",
        "pathloss",
        "barshadow",
        "cal",
        "s2d",
        "x1d",
    ],
)
def test_nirspec_mos_mt_spec2(run_spec2_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test for calwebb_spec2 on a NIRSpec MOS Moving Target exposure."""

    # Run the pipeline and retrieve outputs
    rtdata = run_spec2_pipeline
    output = f"jw01444003001_04102_00001_nrs1_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nirspec_mos_movingtarget/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
@pytest.mark.parametrize("source_id", ["v000000001"])
def test_nirspec_mos_mt_spec3(run_spec3_pipeline, fitsdiff_default_kwargs, suffix, source_id):
    """Regression test for calwebb_spec3 on a NIRSpec MOS Moving Target exposure."""

    # Run the pipeline and retrieve outputs
    rtdata = run_spec3_pipeline
    output = f"jw01444-o003_{source_id}_nirspec_clear-prism_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nirspec_mos_movingtarget/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
