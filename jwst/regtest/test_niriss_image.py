""" Test for the detector1 pipeline using NIRISS image mode, starting with
    an uncal file. The charge_migration and ramp fitting output
    products are saved for comparisons for those two steps.
"""

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_detector1(rtdata_module):
    """Run calwebb_detector1 pipeline on NIRISS imaging data."""
    rtdata = rtdata_module

    rtdata.get_data("niriss/jw01094001002_02107_00001_nis_uncal.fits")

    # Run detector1 pipeline on an _uncal files
    args = ["calwebb_detector1", rtdata.input,
            "--steps.charge_migration.skip=False",
            "--steps.charge_migration.save_results=True",
            "--steps.ramp_fit.save_results=True",
            "--steps.persistence.save_trapsfilled=False",
            ]

    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["charge_migration", "rate", "rateints"])
def test_niriss_image_detector1(run_detector1, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 pipeline performed on NIRISS imaging data.
    """
    _assert_is_same(rtdata_module, fitsdiff_default_kwargs, suffix)


def _assert_is_same(rtdata_module, fitsdiff_default_kwargs, suffix):
    """Assertion helper for the above tests"""
    rtdata = rtdata_module
    rtdata.input = "jw01094001002_02107_00001_nis_uncal.fits"
    output = f"jw01094001002_02107_00001_nis_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_niriss_image/{output}")

    # Set tolerances so the crf, rscd and rateints file comparisons work across
    # architectures
    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 1e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
