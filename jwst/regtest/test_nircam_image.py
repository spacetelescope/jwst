from glob import glob

import pytest
from astropy.io.fits.diff import FITSDiff
from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose

from jwst.stpipe import Step
from jwst import datamodels


@pytest.fixture(scope="module")
def run_detector1pipeline(jail, rtdata_module):
    """Run calwebb_detector1 on NIRCam imaging long data"""
    rtdata = rtdata_module
    rtdata.get_data("nircam/image/jw42424001001_01101_00001_nrca5_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files
    args = ["calwebb_detector1", rtdata.input,
            "--steps.dq_init.save_results=True",
            "--steps.saturation.save_results=True",
            "--steps.superbias.save_results=True",
            "--steps.refpix.save_results=True",
            "--steps.linearity.save_results=True",
            "--steps.dark_current.save_results=True",
            "--steps.jump.save_results=True",
            "--steps.jump.rejection_threshold=50.0",
            ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_image2pipeline(run_detector1pipeline, jail, rtdata_module):
    """Run calwebb_image2 on NIRCam imaging long data"""
    rtdata = rtdata_module
    rtdata.input = "jw42424001001_01101_00001_nrca5_rate.fits"
    args = ["calwebb_image2", rtdata.input,
            "--steps.assign_wcs.save_results=True",
            "--steps.flat_field.save_results=True",
            ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_image3pipeline(run_image2pipeline, rtdata_module, jail):
    """Run calwebb_image3 on NIRCam imaging long data"""
    rtdata = rtdata_module
    # Grab rest of _rate files for the asn and run image2 pipeline on each to
    # produce fresh _cal files for the image3 pipeline.  We won't check these
    # or look at intermediate products, including the resampled i2d
    rate_files = [
        "nircam/image/jw42424001001_01101_00001_nrcb5_rate.fits",
        "nircam/image/jw42424001001_01101_00002_nrca5_rate.fits",
        "nircam/image/jw42424001001_01101_00002_nrcb5_rate.fits",
        "nircam/image/jw42424001001_01101_00003_nrca5_rate.fits",
        "nircam/image/jw42424001001_01101_00003_nrcb5_rate.fits",
    ]
    for rate_file in rate_files:
        rtdata.get_data(rate_file)
        args = ["calwebb_image2", rtdata.input]
        Step.from_cmdline(args)

    # Get the level3 association json file (though not its members) and run
    # image3 pipeline on all _cal files listed in association
    rtdata.get_data("nircam/image/jw42424-o002_20191220t214154_image3_001_asn.json")
    args = ["calwebb_image3", rtdata.input,
            "--steps.tweakreg.save_results=True",
            # "--steps.skymatch.save_results=True",
            "--steps.source_catalog.snr_threshold=20",
            ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["dq_init", "saturation", "superbias",
                                    "refpix", "linearity", "trapsfilled",
                                    "dark_current", "jump", "rate",
                                    "flat_field", "cal", "i2d"])
def test_nircam_image_stages12(run_image2pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 and image2 pipelines performed on NIRCam data."""
    rtdata = rtdata_module
    rtdata.input = "jw42424001001_01101_00001_nrca5_uncal.fits"
    output = "jw42424001001_01101_00001_nrca5_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "i2d":
        fitsdiff_default_kwargs["rtol"] = 5e-5
        fitsdiff_default_kwargs["atol"] = 1e-4

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_image_stage2_wcs(run_image2pipeline, rtdata_module):
    """Test that WCS object works as expected"""
    rtdata = rtdata_module
    rtdata.input = "jw42424001001_01101_00001_nrca5_uncal.fits"
    output = "jw42424001001_01101_00001_nrca5_assign_wcs.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    with datamodels.open(rtdata.output) as model, datamodels.open(rtdata.truth) as model_truth:
        grid = grid_from_bounding_box(model.meta.wcs.bounding_box)

        ra, dec = model.meta.wcs(*grid)
        ra_truth, dec_truth = model_truth.meta.wcs(*grid)

        assert_allclose(ra, ra_truth)
        assert_allclose(dec, dec_truth)


@pytest.mark.bigdata
def test_nircam_image_stage3_tweakreg(run_image3pipeline):
    """Test that tweakreg doesn't attach a catalog and that it updates the wcs"""
    files = glob("*tweakreg.fits")
    for filename in files:
        with datamodels.open(filename) as model:
            # Makes sure the catalog is not attached
            with pytest.raises(AttributeError):
                model.catalog

            # Check that all but the first exposure in the association
            # has a WCS correction applied
            if "jw42424001001_01101_00001" not in model.meta.filename:
                assert "v2v3corr" in model.meta.wcs.available_frames


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["i2d"])
def test_nircam_image_stage3(run_image3pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test that resampled i2d looks good for NIRCam imaging"""
    rtdata = rtdata_module
    rtdata.input = "jw42424-o002_20191220t214154_image3_001_asn.json"
    output = f"jw42424-o002_t001_nircam_clear-f444w_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 2e-4

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_image_stage3_catalog(run_image3pipeline, rtdata_module, diff_astropy_tables):
    rtdata = rtdata_module
    rtdata.input = "jw42424-o002_20191220t214154_image3_001_asn.json"
    output = "jw42424-o002_t001_nircam_clear-f444w_cat.ecsv"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth, rtol=1e-4, atol=1e-5)


@pytest.mark.bigdata
def test_nircam_image_stage3_segm(run_image3pipeline, rtdata_module, fitsdiff_default_kwargs):
    """Test that segmentation map looks good for NIRCam imaging"""
    rtdata = rtdata_module
    rtdata.input = "jw42424-o002_20191220t214154_image3_001_asn.json"
    output = "jw42424-o002_t001_nircam_clear-f444w_segm.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.fixture()
def run_image3_closedfile(rtdata, jail):
    """Run calwebb_image3 on NIRCam imaging with data that had a closed file issue."""
    rtdata.get_asn("nircam/image/fail_short_image3_asn.json")

    args = ["calwebb_image3", rtdata.input]
    Step.from_cmdline(args)


@pytest.mark.bigdata
def test_image3_closedfile(run_image3_closedfile, rtdata, fitsdiff_default_kwargs):
    """Ensure production of Image3Pipeline output with data having closed file issues"""
    rtdata.output = 'jw00617-o082_t001_nircam_clear-f090w-sub320_i2d.fits'
    rtdata.get_truth('truth/test_nircam_image/jw00617-o082_t001_nircam_clear-f090w-sub320_i2d.fits')

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_frame_averaged_darks(rtdata, fitsdiff_default_kwargs):
    """Test optional frame-averaged darks output from DarkCurrentStep"""
    rtdata.get_data("nircam/image/jw00312007001_02102_00001_nrcblong_ramp.fits")

    args = ["jwst.dark_current.DarkCurrentStep", rtdata.input,
            "--dark_output='frame_averaged_darks.fits'",
            ]
    Step.from_cmdline(args)
    rtdata.output = "frame_averaged_darks.fits"

    rtdata.get_truth("truth/test_nircam_image/frame_averaged_darks.fits")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
