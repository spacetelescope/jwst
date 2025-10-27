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
    spec2_asns = [
        "miri/wfss/jw03224-test_spec2_00001_asn.json",
    ]

    foffset_file = "jwst_miri_filteroffset_0008.asdf"
    photom_file = "jwst_miri_photom_0217.fits"
    extract1d_file = "jwst_miri_extract1d_0006.json"
    distortion_file = "jwst_miri_distortion_0047.asdf"

    wave_file = "miri_wfss_wavelengthrange.asdf"
    apcorr_file = "jwst_miri_apcorr_0015.fits"
    rtdata.get_data(f"miri/wfss/{foffset_file}")
    rtdata.get_data(f"miri/wfss/{photom_file}")
    rtdata.get_data(f"miri/wfss/{extract1d_file}")
    rtdata.get_data(f"miri/wfss/{distortion_file}")
    rtdata.get_data(f"miri/wfss/{wave_file}")
    rtdata.get_data(f"miri/wfss/{apcorr_file}")

    rtdata.get_asn(spec2_asns[0])
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        f"--steps.assign_wcs.override_filteroffset={foffset_file}",
        f"--steps.assign_wcs.override_distortion={distortion_file}",
        "--steps.assign_wcs.override_regions=N/A",
        f"--steps.assign_wcs.override_wavelengthrange={wave_file}",
        f"--steps.photom.override_photom={photom_file}",
        f"--steps.extract_1d.override_extract1d={extract1d_file}",
        f"--steps.extract_1d.override_apcorr={apcorr_file}",
        "--steps.bkg_subtract.skip=true",
        "--steps.flat_field.save_results=true",
        "--steps.extract_2d.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.photom.save_results=true",
        "--steps.extract_1d.save_results=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.mark.parametrize(
    "suffix",
    ["assign_wcs", "cal", "extract_2d", "flat_field", "photom", "srctype", "x1d"],
)
def test_miri_wfss_spec2(run_miri_wfss_spec2, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test for calwebb_spec2 applied to MIRI WFSS data"""
    rtdata = rtdata_module
    rtdata.input = "jw03224024001_03103_00001_mirimage_rate.fits"
    output = "jw03224024001_03103_00001_mirimage_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_miri_wfss/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
