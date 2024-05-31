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
    rtdata.get_asn("nirspec/mos/jw01345-o066_20230831t181155_spec3_00002_asn.json")

    # Run the calwebb_spec3 pipeline on the association
    args = ["calwebb_spec3", rtdata.input]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
@pytest.mark.parametrize("source_id", ["b000000030", "b000000031",
                                       "s000004385", "s000007380",
                                       "v000000048", "v000000049",
                                       "v000000053", "v000000056"])
def test_nirspec_mos_spec3(run_pipeline, suffix, source_id, fitsdiff_default_kwargs):
    """Check results of calwebb_spec3"""
    rtdata = run_pipeline

    output = f"jw01345-o066_{source_id}_nirspec_f170lp-g235m_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_mos_spec3/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "s2d":
        fitsdiff_default_kwargs["rtol"] = 1e-4
        fitsdiff_default_kwargs["atol"] = 1e-5

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    if suffix == "s2d":
        # Compare the calculated wavelengths
        tolerance = 1e-03
        dmt = datamodels.open(rtdata.truth)
        dmr = datamodels.open(rtdata.output)
        if isinstance(dmt, datamodels.MultiSlitModel):
            names = [s.name for s in dmt.slits]
            for name in names:
                st_idx = [(s.wcs, s.wavelength) for s in dmt.slits if s.name==name]
                w = dmt.slits[st_idx].meta.wcs
                x, y = wcstools.grid_from_bounding_box(w.bounding_box, step=(1, 1), center=True)
                _, _, wave = w(x, y)
                sr_idx = [(s.wcs, s.wavelength) for s in dmr.slits if s.name==name]
                wlr = dmr.slits[sr_idx].wavelength
                assert np.all(np.isclose(wave, wlr, atol=tolerance))
        else:
            w = dmt.meta.wcs
            x, y = wcstools.grid_from_bounding_box(w.bounding_box, step=(1, 1), center=True)
            _, _, wave = w(x, y)
            wlr = dmr.wavelength
            assert np.all(np.isclose(wave, wlr, atol=tolerance))
