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

    photom_file = "drlphotom.fits"
    flat_file = "jwst_miri_flat_0850.fits"
    filteroffset_file = "jwst_miri_filteroffset_0019.asdf"
    specwcs_file = "MIRI_WFSS_specwcs_20260316.asdf"
    rtdata.get_data(f"miri/wfss/{photom_file}")
    rtdata.get_data(f"miri/wfss/{flat_file}")
    rtdata.get_data(f"miri/wfss/{filteroffset_file}")
    rtdata.get_data(f"miri/wfss/{specwcs_file}")

    rtdata.get_asn("miri/wfss/jw09505-o001_spec2_00001_asn.json")
    args = [
        "calwebb_spec2",
        rtdata.input,
        f"--steps.photom.override_photom={photom_file}",
        f"--steps.flat_field.override_flat={flat_file}",
        f"--steps.assign_wcs.override_specwcs={specwcs_file}",
        f"--steps.assign_wcs.override_filteroffset={filteroffset_file}",
        "--steps.bkg_subtract.skip=false",
        "--steps.flat_field.skip=false",
        "--steps.assign_wcs.save_results=true",
        "--steps.extract_2d.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.photom.save_results=true",
        "--steps.flat_field.save_results=true",
        "--steps.bkg_subtract.save_results=true",
        "--steps.extract_1d.save_results=true",
    ]

    Step.from_cmdline(args)


@pytest.mark.parametrize(
    "suffix",
    ["assign_wcs", "cal", "extract_2d", "photom", "srctype", "x1d", "bsub", "flat_field"],
)
def test_miri_wfss_spec2(run_miri_wfss_spec2, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test for calwebb_spec2 applied to MIRI WFSS data"""
    rtdata = rtdata_module
    rtdata.input = "jw09505001001_02102_00004_mirimage_rate.fits"
    output = "jw09505001001_02102_00004_mirimage_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_miri_wfss/{output}")

    # Ignore the custom extract1d file because it contains a full path.
    fitsdiff_default_kwargs["ignore_keywords"].append("R_PHOTOM")
    fitsdiff_default_kwargs["ignore_keywords"].append("R_FLAT")
    fitsdiff_default_kwargs["ignore_keywords"].append("R_SPCWCS")
    fitsdiff_default_kwargs["ignore_keywords"].append("R_FILOFF")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
