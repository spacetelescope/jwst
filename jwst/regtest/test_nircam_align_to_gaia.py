import pytest
from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose

from stdatamodels.jwst import datamodels
from jwst.stpipe import Step
from jwst.tweakreg import TweakRegStep


@pytest.fixture(scope="module")
def run_image3pipeline(rtdata_module):
    """Run calwebb_image3 on NIRCam imaging and align to gaia"""

    rtdata = rtdata_module
    rtdata.get_asn("nircam/image/level3_F277W_3img_asn.json")
    args = [
        "calwebb_image3",
        rtdata.input,
        "--steps.tweakreg.abs_refcat=GAIADR2",
        "--steps.tweakreg.save_results=True",
        "--steps.tweakreg.output_use_model=True",
    ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "root", ["jw01069002001_01101_00001", "jw01069002003_01101_00009", "jw01069002004_01101_00013"]
)
def test_tweakreg_with_gaia(run_image3pipeline, rtdata_module, root):
    """Test that WCS object works as expected"""
    rtdata = rtdata_module
    rtdata.input = root + "_nrca5_cal.fits"
    output_crf = root + "_nrca5_a3001_crf.fits"
    rtdata.output = output_crf
    rtdata.get_truth("truth/test_nircam_align_to_gaia/" + output_crf)

    with datamodels.open(rtdata.output) as model, datamodels.open(rtdata.truth) as model_truth:
        grid = grid_from_bounding_box(model.meta.wcs.bounding_box)

        ra, dec = model.meta.wcs(*grid)
        ra_truth, dec_truth = model_truth.meta.wcs(*grid)

        assert_allclose(ra, ra_truth)
        assert_allclose(dec, dec_truth)


@pytest.mark.bigdata
def test_sourcecat_as_abs_refcat(run_image3pipeline, rtdata_module):
    """
    Test that source catalog output is compatible with tweakreg.

    Some workflows require using source catalog step output as the
    absolute reference catalog for tweakreg. This test ensures that
    this is possible. This is not nircam-specific, but is included here
    to avoid adding an additional test run through image3.
    """
    rtdata = rtdata_module
    rtdata.output = "LMC_F277W_modA_dither_mosaic_cat.ecsv"  # output from source catalog step
    rtdata.get_asn("nircam/image/level3_F277W_3img_asn.json")

    TweakRegStep.call(rtdata.input, abs_refcat=rtdata.output)
