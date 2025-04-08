import pytest

from astropy.io.fits.diff import FITSDiff

from jwst.regtest.regtestdata import RELAX_TOL
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run the calwebb_spec2 pipeline on a single NIRSpec WFSS exposure."""

    rtdata = rtdata_module

    # Get the input data; load individual data files first, load ASN file last
    rtdata.get_data("nircam/wfss/jw01076-o101_t002_nircam_clear-f356w_cat.ecsv")
    rtdata.get_data("nircam/wfss/jw01076101001_02101_00003_nrcalong_rate.fits")
    rtdata.get_data("nircam/wfss/jw01076-o101_20220403t120233_spec2_002_asn.json")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = ["calwebb_spec2", rtdata.input,
            "--steps.assign_wcs.save_results=true",
            "--steps.bkg_subtract.skip=False",
            "--steps.bkg_subtract.wfss_mmag_extract=20.0",
            "--steps.bkg_subtract.save_results=true",
            "--steps.extract_2d.save_results=true",
            "--steps.extract_2d.wfss_mmag_extract=19.0",
            "--steps.extract_2d.wfss_nbright=20",
            "--steps.srctype.save_results=true",
            "--steps.flat_field.save_results=true"]
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
    if RELAX_TOL and suffix == "bsub":
        # 1671 different pixels found (0.04% different).
        # Maximum relative difference: 0.00041919932118617
        # Maximum absolute difference: 1.1920928955078125e-07
        fitsdiff_default_kwargs["rtol"] = 5e-4
        fitsdiff_default_kwargs["atol"] = 2e-7
    elif RELAX_TOL and suffix == "extract_2d":
        # 16 different pixels found (0.03% different).
        # Maximum relative difference: 0.0020618431735783815
        # Maximum absolute difference: 1.2549571692943573e-07
        fitsdiff_default_kwargs["rtol"] = 0.003
        fitsdiff_default_kwargs["atol"] = 2e-7
    elif RELAX_TOL and suffix == "flat_field":
        # 2614 different pixels found (0.06% different).
        # Maximum relative difference: 0.0283019058406353
        # Maximum absolute difference: 1.2980308383703232e-07
        fitsdiff_default_kwargs["rtol"] = 0.03
        fitsdiff_default_kwargs["atol"] = 2e-7
    elif RELAX_TOL and suffix == "srctype":
        # 48 different pixels found (0.05% different).
        # Maximum relative difference: 0.00050398672465235
        # Maximum absolute difference: 1.310836523771286e-07
        fitsdiff_default_kwargs["rtol"] = 6e-4
        fitsdiff_default_kwargs["atol"] = 2e-7
    elif RELAX_TOL and suffix == "cal":
        # Each extension has different failures so using max from combo tol.
        fitsdiff_default_kwargs["rtol"] = 0.2
        fitsdiff_default_kwargs["atol"] = 8
    elif RELAX_TOL and suffix == "x1d":
        # Each extension has different failures so using max from combo tol.
        fitsdiff_default_kwargs["rtol"] = 0.01
        fitsdiff_default_kwargs["atol"] = 0.2
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
