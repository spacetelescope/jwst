"""Regression tests for MIRI WFSS mode"""

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_miri_wfss_spec2(rtdata_module, resource_tracker):
    """Run the calwebb_spec2 pipeline on MIRI WFSS exposures"""
    rtdata = rtdata_module
    # These are the WFSS exposures we'll be processing

    photom_file = "jwst_miri_photom_WFSS_20260220_v2.fits"
    bkg_file = "MIRI_WFSS_bkg_February2026.fits"
    flat_file = "jwst_miri_flat_WFSS_20260220.fits"
    rtdata.get_data(f"miri/wfss/{photom_file}")
    rtdata.get_data(f"miri/wfss/{bkg_file}")
    rtdata.get_data(f"miri/wfss/{flat_file}")

    rtdata.get_asn("miri/wfss/jw03224-miri_wfss_spec2_00001_asn.json")
    args = [
        "calwebb_spec2",
        rtdata.input,
        f"--steps.photom.override_photom={photom_file}",
        f"--steps.flat.override_flat={flat_file}",
        f"--steps.bkg_subtract.override_bkg_sutract={bkg_file}",
        "--steps.bkg_subtract.skip=false",
        "--steps.flat_field.skip=false",
        "--steps.assign_wcs.save_results=true",
        "--steps.extract_2d.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.photom.save_results=true",
        "--steps.flat.save_results=true",
        "--steps.bkg_subtract.save_results=true",
        "--steps.extract_1d.save_results=true",
    ]

    Step.from_cmdline(args)


@pytest.mark.parametrize(
    "suffix",
    ["assign_wcs", "cal", "extract_2d", "photom", "srctype", "x1d", "bsub", "flat"],
)
def test_miri_wfss_spec2(run_miri_wfss_spec2, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test for calwebb_spec2 applied to MIRI WFSS data"""
    rtdata = rtdata_module
    rtdata.input = "jw03224024001_03103_00001_mirimage_rate.fits"
    output = "jw03224_miri_wfss_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_miri_wfss/{output}")

    # Ignore the custom extract1d file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_PHOTOM")
    fitsdiff_default_kwargs["ignore_keywords"].append("R_SPCWCS")
    fitsdiff_default_kwargs["ignore_keywords"].append("R_DISTOR")
    fitsdiff_default_kwargs["ignore_keywords"].append("R_WAVRAN")
    fitsdiff_default_kwargs["ignore_keywords"].append("R_FILOFF")
    fitsdiff_default_kwargs["ignore_keywords"].append("R_EXTR1D")
    fitsdiff_default_kwargs["ignore_keywords"].append("R_REGION")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
