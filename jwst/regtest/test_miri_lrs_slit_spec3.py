""" Test of the spec3 pipeline using MIRI LRS fixed-slit exposures.
    This takes an association and generates the level 3 products."""
import pytest
import numpy as np
from gwcs import wcstools

import asdf
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step
from jwst import datamodels


@pytest.fixture(scope='module')
def run_pipeline(rtdata_module, request):
    """
    Run the calwebb_spec3 pipeline on an ASN of nodded MIRI LRS
    fixed-slit exposures.
    """
    rtdata = rtdata_module

    # Get the spec3 ASN and its members
    rtdata.get_asn("miri/lrs/jw01530-o005_20221202t204827_spec3_00001_asn.json")
    args = [
        "calwebb_spec3",
        rtdata.input
    ]

    # Run the calwebb_spec3 pipeline; save results from intermediate steps
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["s2d", "x1d"])
def test_miri_lrs_slit_spec3(run_pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec3 pipeline on MIRI
       LRS fixed-slit data using along-slit-nod pattern for
       background subtraction."""

    # Run the pipeline and retrieve outputs
    rtdata = rtdata_module

    output = f"jw01530-o005_t004_miri_p750l_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(f"truth/test_miri_lrs_slit_spec3/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    if output == "s2d":
        # Compare the calculated wavelengths
        tolerance = 1e-03
        dmt = datamodels.open(rtdata.truth)
        dmr = datamodels.open(rtdata.output)
        if isinstance(dmt, datamodels.MultiSlitModel):
            names = [s.name for s in dmt.slits]
            for name in names:
                st_idx = [(s.wcs, s.wavelength) for s in dmt.slits if s.name == name]
                w = dmt.slits[st_idx].meta.wcs
                x, y = wcstools.grid_from_bounding_box(w.bounding_box, step=(1, 1), center=True)
                _, _, wave = w(x, y)
                sr_idx = [(s.wcs, s.wavelength) for s in dmr.slits if s.name == name]
                wlr = dmr.slits[sr_idx].wavelength
                assert np.all(np.isclose(wave, wlr, atol=tolerance))
        else:
            w = dmt.meta.wcs
            x, y = wcstools.grid_from_bounding_box(w.bounding_box, step=(1, 1), center=True)
            _, _, wave = w(x, y)
            wlr = dmr.wavelength
            assert np.all(np.isclose(wave, wlr, atol=tolerance))
