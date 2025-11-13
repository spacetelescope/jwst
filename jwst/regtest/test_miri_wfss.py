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

    foffset_file = "jwst_miri_filteroffset_0008.asdf"
    photom_file = "jwst_miri_photom_0217.fits"
    extract1d_file = "jwst_miri_extract1d_0006.json"
    distortion_file = "jwst_miri_distortion_0047.asdf"
    specwcs_file = "MIRI_WFSS_specwcs_20250911.asdf"
    wave_file = "miri_wfss_wavelengthrange.asdf"
    apcorr_file = "jwst_miri_apcorr_0015.fits"
    rtdata.get_data(f"miri/wfss/{foffset_file}")
    rtdata.get_data(f"miri/wfss/{photom_file}")
    rtdata.get_data(f"miri/wfss/{extract1d_file}")
    rtdata.get_data(f"miri/wfss/{distortion_file}")
    rtdata.get_data(f"miri/wfss/{wave_file}")
    rtdata.get_data(f"miri/wfss/{apcorr_file}")
    rtdata.get_data(f"miri/wfss/{specwcs_file}")

    rtdata.get_asn("miri/wfss/jw03224-miri_wfss_spec2_00001_asn.json")
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.assign_wcs.skip=false",
        f"--steps.assign_wcs.override_filteroffset={foffset_file}",
        f"--steps.assign_wcs.override_distortion={distortion_file}",
        f"--steps.assign_wcs.override_specwcs={specwcs_file}",
        "--steps.assign_wcs.override_regions=N/A",
        f"--steps.assign_wcs.override_wavelengthrange={wave_file}",
        f"--steps.extract_2d.override_wavelengthrange={wave_file}",
        f"--steps.photom.override_photom={photom_file}",
        f"--steps.extract_1d.override_extract1d={extract1d_file}",
        f"--steps.extract_1d.override_apcorr={apcorr_file}",
        "--steps.bkg_subtract.skip=true",
        "--steps.flat_field.skip=true",
        "--steps.extract_2d.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.photom.save_results=true",
        "--steps.extract_1d.save_results=true",
    ]

    Step.from_cmdline(args)


@pytest.mark.parametrize(
    "suffix",
    ["assign_wcs", "cal", "extract_2d", "photom", "srctype", "x1d"],
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
