import pytest
from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose

from jwst.stpipe import Step
from jwst import datamodels


@pytest.fixture(scope="module")
def run_image3pipeline(rtdata_module, jail):
    ''' Run calwebb_image3 on NIRCam imaging and align to gaia '''

    rtdata = rtdata_module
    rtdata.get_asn("nircam/image/level3_F277W_3img_asn.json")
    args = ["calwebb_image3", rtdata.input,
            "--steps.tweakreg.abs_refcat=GAIADR2",
            "--steps.tweakreg.save_results=True",
            "--steps.tweakreg.output_use_model=True"
            ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("root", ["jw01069002001_01101_00001",
                                  "jw01069002003_01101_00009",
                                  "jw01069002004_01101_00013"])
def test_tweakreg_with_gaia(run_image3pipeline, rtdata_module, root):
    """ Test that WCS object works as expected """
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
