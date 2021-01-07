import pytest
from astropy.io.fits.diff import FITSDiff
from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step
from jwst import datamodels


@pytest.fixture(scope="module")
def run_detector1pipeline(jail, rtdata_module):
    """Run calwebb_detector1 on NIRCam imaging long data"""
    rtdata = rtdata_module
    rtdata.get_data("nircam/image/jw42424001001_01101_00001_nrca5_uncal.fits")

    collect_pipeline_cfgs("config")

    # Run detector1 pipeline only on one of the _uncal files
    args = ["config/calwebb_detector1.cfg", rtdata.input,
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
    args = ["config/calwebb_image2.cfg", rtdata.input,
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
        args = ["config/calwebb_image2.cfg", rtdata.input,
            "--steps.resample.skip=True"]
        Step.from_cmdline(args)

    # Get the level3 assocation json file (though not its members) and run
    # image3 pipeline on all _cal files listed in association
    rtdata.get_data("nircam/image/jw42424-o002_20191220t214154_image3_001_asn.json")
    args = ["config/calwebb_image3.cfg", rtdata.input,
        # Comment out following lines, as the dataset is currently broken
        # "--steps.tweakreg.save_results=True",
        # "--steps.skymatch.save_results=True",
        "--steps.source_catalog.snr_threshold=20",
        ]
    Step.from_cmdline(args)


@pytest.fixture()
def run_image3_closedfile(rtdata, jail):
    """Run calwebb_image3 on NIRCam imaging with data that had a closed file issue."""

    rtdata.get_asn("nircam/image/fail_short_image3_asn.json")

    collect_pipeline_cfgs("config")

    args = ["config/calwebb_image3.cfg", rtdata.input]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["dq_init", "saturation", "superbias",
    "refpix", "linearity", "trapsfilled", "dark_current", "jump", "rate",
    "flat_field", "cal", "i2d"])
def test_nircam_image_stages12(run_image2pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 and image2 pipelines performed on NIRCam data."""
    rtdata = rtdata_module
    rtdata.input = "jw42424001001_01101_00001_nrca5_uncal.fits"
    output = "jw42424001001_01101_00001_nrca5_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nircam_image_stages/" + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_image_stage2_wcs(run_image2pipeline, rtdata_module):
    """Test that WCS object works as expected"""
    rtdata = rtdata_module
    rtdata.input = "jw42424001001_01101_00001_nrca5_uncal.fits"
    output = "jw42424001001_01101_00001_nrca5_assign_wcs.fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nircam_image_stages/" + output)

    with datamodels.open(rtdata.output) as model, datamodels.open(rtdata.truth) as model_truth:
        grid = grid_from_bounding_box(model.meta.wcs.bounding_box)

        ra, dec = model.meta.wcs(*grid)
        ra_truth, dec_truth = model_truth.meta.wcs(*grid)

        assert_allclose(ra, ra_truth)
        assert_allclose(dec, dec_truth)


@pytest.mark.bigdata
def test_nircam_image_stage3_i2d(run_image3pipeline, rtdata_module, fitsdiff_default_kwargs):
    """Test that resampled i2d looks good for NIRCam imaging"""
    rtdata = rtdata_module
    rtdata.input = "jw42424-o002_20191220t214154_image3_001_asn.json"
    rtdata.output = "jw42424-o002_t001_nircam_clear-f444w_i2d.fits"
    rtdata.get_truth("truth/test_nircam_image_stages/jw42424-o002_t001_nircam_clear-f444w_i2d.fits")

    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_image_stage3_catalog(run_image3pipeline, rtdata_module, diff_astropy_tables):
    rtdata = rtdata_module
    rtdata.input = "jw42424-o002_20191220t214154_image3_001_asn.json"
    rtdata.output = "jw42424-o002_t001_nircam_clear-f444w_cat.ecsv"
    rtdata.get_truth("truth/test_nircam_image_stages/jw42424-o002_t001_nircam_clear-f444w_cat.ecsv")

    assert diff_astropy_tables(rtdata.output, rtdata.truth, rtol=1e-4, atol=1e-5)


@pytest.mark.bigdata
def test_image3_closedfile(run_image3_closedfile, rtdata, fitsdiff_default_kwargs):
    """Ensure production of Image3Pipeline output with data having closed file issues"""
    rtdata.output = 'jw00617-o082_t001_nircam_clear-f090w-sub320_i2d.fits'
    rtdata.get_truth('truth/test_nircam_image/jw00617-o082_t001_nircam_clear-f090w-sub320_i2d.fits')

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
