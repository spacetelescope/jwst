"""Test of the spec3 pipeline using MIRI LRS fixed-slit exposures.
This takes an association and generates the level 3 products."""

import pytest
import numpy as np


from gwcs import wcstools
import asdf
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.stpipe import Step
from jwst import datamodels

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(
    scope="module", params=["default_wcs", "user_wcs", "user_wcs+shape", "user_wcs+shape1"]
)
def run_pipeline(rtdata_module, request, resource_tracker):
    """
    Run the calwebb_spec3 pipeline on an ASN of nodded MIRI LRS
    fixed-slit exposures using different options for the WCS and output
    image shape for the resample step.

    The first iteration ("default_wcs") creates an output WCS for the combined
    product based on the built-in WCS's in the inputs.

    The second iteration ("user_wcs") creates a user-supplied WCS input file
    using the WCS from the default product, just to prove that you get the
    identical result when using a user-supplied WCS.

    The third iteration ("user_wcs+shape") uses the same user-supplied WCS and also
    specifies the output 2D file shape, using a shape that is identical to the
    default. Hence all of the first 3 iterations should produce identical results
    and therefore are compared to a single set of truth files.

    The fourth iteration ("user_wcs+shape1") uses the same user-specified WCS and
    specifies an output 2D file shape that is 1 pixel larger than the default in
    both axes. Hence the resulting s2d product needs a separate (larger) truth file.
    Meanwhile, the x1d product from this iteration should still be identical to
    the first 3, because the extra row and column of the 2D data file are ignored
    during extraction.
    """
    rtdata = rtdata_module

    # Get the spec3 ASN and its members
    rtdata.get_asn("miri/lrs/jw01530-o005_20221202t204827_spec3_00001_asn.json")
    root_file = "jw01530-o005_t004_miri_p750l_"

    args = ["calwebb_spec3", rtdata.input]

    rtdata.custom_wcs_mode = request.param

    if request.param != "default_wcs":
        # Get the s2d product that was just created using "default_wcs"
        default_s2d = root_file + "s2d.fits"
        dm = datamodels.open(default_s2d)
        # Create a user-supplied WCS file that is identical to the default WCS
        af = asdf.AsdfFile({"wcs": dm.meta.wcs})
        wcs_file = default_s2d[:-8] + "wcs.asdf"
        af.write_to(wcs_file)
        args.append(f"--steps.resample_spec.output_wcs={wcs_file}")

    if request.param == "user_wcs+shape":
        output_shape = ",".join(map(str, dm.data.shape[::-1]))
        args.append(f"--steps.resample_spec.output_shape={output_shape}")

    elif request.param == "user_wcs+shape1":
        output_shape = ",".join(map(str, (d + 1 for d in dm.data.shape[::-1])))
        args.append(f"--steps.resample_spec.output_shape={output_shape}")
        output_file = root_file + "shape1.fits"
        args.append(f"--steps.resample_spec.output_file={output_file}")

    # Run the calwebb_spec3 pipeline; save results from intermediate steps
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_spec3(log_tracked_resources, run_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["s2d", "x1d"])
def test_miri_lrs_slit_spec3(run_pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec3 pipeline on MIRI
    LRS fixed-slit data using along-slit-nod pattern for
    background subtraction."""

    # Run the pipeline and retrieve outputs
    rtdata = rtdata_module

    if rtdata.custom_wcs_mode == "user_wcs+shape1" and suffix == "s2d":
        output = f"jw01530-o005_t004_miri_p750l_shape1_{suffix}.fits"
    else:
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
