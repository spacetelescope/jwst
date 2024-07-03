from astropy.io.fits.diff import FITSDiff
import pytest
import numpy as np
from gwcs import wcstools

from jwst.stpipe import Step
from stdatamodels.jwst import datamodels


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module):
    """
    Run the calwebb_spec3 pipeline on a NIRSpec FS moving target.
    """
    rtdata = rtdata_module

    # Get the ASN file and input exposures
    rtdata.get_asn('nirspec/fs/jw01245-o002_20240701t053319_spec3_00001_asn.json')

    # Run the calwebb_spec3 pipeline; save results from intermediate steps
    args = ["calwebb_spec3", rtdata.input,
            "--steps.outlier_detection.save_results=true",
            "--steps.resample_spec.save_results=true",
            "--steps.extract_1d.save_results=true"]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
def test_nirspec_fs_spec3_moving_target(
        run_pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test spec3 pipeline on a NIRSpec FS moving target."""
    rtdata = rtdata_module

    output = f"jw01245-o002_s000000001_nirspec_clear-prism-s200a1-subs200a1_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_fs_spec3_moving_target/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "s2d":
        fitsdiff_default_kwargs["rtol"] = 1e-2
        fitsdiff_default_kwargs["atol"] = 2e-4

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
