""" Test of the spec3 pipeline using MIRI LRS fixed-slit exposures.
    This takes an association and generates the level 3 products."""
import pytest

import asdf
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step
from jwst import datamodels


@pytest.fixture(
    scope="module",
    params=["default_wcs", "user_wcs", "user_wcs+shape", "user_wcs+shape1"]
)
def run_pipeline(jail, rtdata_module, request):
    """
    Run the calwebb_spec3 pipeline on an ASN of nodded MIRI LRS
    fixed-slit exposures using different options for the WCS and output
    image shape for the resample step.
    """
    rtdata = rtdata_module

    # Get the spec3 ASN and its members
    rtdata.get_asn("miri/lrs/jw01530-o005_20221202t204827_spec3_00001_asn.json")
    root_file = "jw01530-o005_t004_miri_p750l_"
    rtdata.get_truth("truth/test_miri_lrs_slit_spec3/jw01530-o005_t004_miri_p750l_s2d.fits")

    args = [
        "calwebb_spec3",
        rtdata.input
    ]

    rtdata.custom_wcs_mode = request.param

    if request.param != "default_wcs":
        dm = datamodels.open(rtdata.truth)
        af = asdf.AsdfFile({"wcs": dm.meta.wcs})
        wcs_file = rtdata.truth[:-8] + 'wcs.asdf'
        af.write_to(wcs_file)
        args.append(f"--steps.resample_spec.output_wcs={wcs_file}")

    if request.param == "user_wcs+shape":
        output_shape = ','.join(map(str, dm.data.shape[::-1]))
        args.append(f"--steps.resample_spec.output_shape={output_shape}")

    elif request.param == "user_wcs+shape1":
        output_shape = ','.join(map(str, (d + 1 for d in dm.data.shape[::-1])))
        args.append(f"--steps.resample_spec.output_shape={output_shape}")
        output_file = root_file + 'shape1.fits'
        args.append(f"--steps.resample_spec.output_file={output_file}")

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

    if rtdata.custom_wcs_mode == 'user_wcs+shape1' and suffix == "s2d":
        output = f"jw01530-o005_t004_miri_p750l_shape1_{suffix}.fits"
    else:
        output = f"jw01530-o005_t004_miri_p750l_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(f"truth/test_miri_lrs_slit_spec3/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
