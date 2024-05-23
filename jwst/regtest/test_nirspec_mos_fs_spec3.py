import pytest
from astropy.io.fits.diff import FITSDiff
import numpy as np
from gwcs import wcstools

from jwst.stpipe import Step
from stdatamodels.jwst import datamodels


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """Run calwebb_spec3 on NIRSpec MOS data."""
    rtdata = rtdata_module
    rtdata.get_asn("nirspec/mos/jw02674-o004_20240305t054741_spec3_00001_asn.json")

    # Run the calwebb_spec3 pipeline on the association
    args = ["calwebb_spec3", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
@pytest.mark.parametrize("source_id", ["s01354", "s12105", "s34946",
                                       "s34949", "s34950", "s34951", "s34952",
                                       "s34953", "s34954", "s34955"])
def test_nirspec_mos_fs_spec3(run_pipeline, suffix, source_id, fitsdiff_default_kwargs):
    """Check results of calwebb_spec3"""
    rtdata = run_pipeline

    output = f"jw02674-o004_{source_id}_nirspec_f290lp-g395m_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_mos_fs_spec3/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "s2d":
        fitsdiff_default_kwargs["rtol"] = 1e-4
        fitsdiff_default_kwargs["atol"] = 1e-5

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
