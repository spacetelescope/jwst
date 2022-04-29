import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run the calwebb_spec2 pipeline on a single NIRSpec WFSS exposure."""

    rtdata = rtdata_module

    # Get the custom pars-spec2pipeline ref file needed to run the pipeline
    rtdata.get_data("nircam/wfss/jwst_nircam_pars-spec2pipeline_custom.asdf")

    # Get the input data; load individual data files first, load ASN file last
    rtdata.get_data("nircam/wfss/jw01076-o101_t002_nircam_clear-f356w_cat.ecsv")
    rtdata.get_data("nircam/wfss/jw01076101001_02101_00003_nrcalong_rate.fits")
    rtdata.get_data("nircam/wfss/jw01076-o101_20220403t120233_spec2_002_asn.json")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    # use "disable-crds-steppars" to exclude pars-spec2pipeline from CRDS
    args = ["jwst_nircam_pars-spec2pipeline_custom.asdf", rtdata.input,
            "--steps.assign_wcs.save_results=true",
            "--steps.bkg_subtract.save_results=true",
            "--steps.extract_2d.save_results=true",
            "--steps.srctype.save_results=true",
            "--steps.flat_field.save_results=true",
            "--disable-crds-steppars"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", [
    "assign_wcs", "bsub", "extract_2d", "flat_field", "srctype",
    "cal", "x1d"])
def test_nircam_wfss_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline on a
       NIRCam WFSS exposure."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = f"jw01076101001_02101_00003_nrcalong_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth("truth/test_nircam_wfss_spec2/" + output)

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
