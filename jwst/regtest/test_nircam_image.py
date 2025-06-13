from glob import glob

import pytest
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose

from stdatamodels.jwst import datamodels

from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_detector1pipeline(rtdata_module):
    """Run calwebb_detector1 on NIRCam imaging long data"""
    rtdata = rtdata_module
    rtdata.get_data("nircam/image/jw01538046001_03105_00001_nrcalong_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files
    args = [
        "calwebb_detector1",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.superbias.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.refpix.refpix_algorithm=median",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.jump.save_results=True",
        "--steps.jump.rejection_threshold=50.0",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_detector1pipeline_with_sirs(rtdata_module, resource_tracker):
    """Run calwebb_detector1 on NIRCam imaging data using SIRS.

    SIRS is the convolution kernel algorithm - Simple Improved Reference Subtraction.
    """
    rtdata = rtdata_module
    rtdata.get_data("nircam/image/jw01345001001_10201_00001_nrca3_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files
    args = [
        "calwebb_detector1",
        rtdata.input,
        "--output_file=jw01345001001_10201_00001_nrca3_sirs",
        "--steps.refpix.refpix_algorithm=sirs",
        "--steps.refpix.save_results=True",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_detector1_with_clean_flicker_noise(rtdata_module):
    """Run detector1 pipeline on NIRCam imaging data with noise cleaning."""
    rtdata_module.get_data("nircam/image/jw01538046001_03105_00001_nrcalong_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files.
    # Run optional clean_flicker_noise step, saving extra outputs
    args = [
        "jwst.pipeline.Detector1Pipeline",
        rtdata_module.input,
        "--output_file=jw01538046001_03105_00001_nrcalong_cfn",
        "--save_calibrated_ramp=True",
        "--steps.clean_flicker_noise.skip=False",
        "--steps.clean_flicker_noise.save_results=True",
        "--steps.clean_flicker_noise.save_mask=True",
        "--steps.clean_flicker_noise.save_background=True",
        "--steps.clean_flicker_noise.save_noise=True",
    ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_detector1_with_likelihood(rtdata_module, resource_tracker):
    """
    Run calwebb_detector1 on NIRCam imaging data using likelihood ramp fitting.

    Input data has a separate zero frame attached.
    """
    rtdata = rtdata_module
    rtdata.get_data("nircam/image/jw01345001001_10201_00001_nrca3_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files
    args = [
        "calwebb_detector1",
        rtdata.input,
        "--output_file=jw01345001001_10201_00001_nrca3_likely",
        "--steps.ramp_fit.algorithm=LIKELY",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_image2pipeline(run_detector1pipeline, rtdata_module, resource_tracker):
    """Run calwebb_image2 on NIRCam imaging long data"""
    rtdata = rtdata_module
    rtdata.input = "jw01538046001_03105_00001_nrcalong_rate.fits"
    args = [
        "calwebb_image2",
        rtdata.input,
        "--steps.assign_wcs.save_results=True",
        "--steps.flat_field.save_results=True",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_image3pipeline(run_image2pipeline, rtdata_module, resource_tracker):
    """Run calwebb_image3 on NIRCam imaging long data"""
    rtdata = rtdata_module
    # Grab rest of _rate files for the asn and run image2 pipeline on each to
    # produce fresh _cal files for the image3 pipeline.  We won't check these
    # or look at intermediate products, including the resampled i2d
    rate_files = [
        "nircam/image/jw01538046001_03105_00001_nrcblong_rate.fits",
        "nircam/image/jw01538046001_03105_00002_nrcalong_rate.fits",
        "nircam/image/jw01538046001_03105_00002_nrcblong_rate.fits",
        "nircam/image/jw01538046001_0310f_00001_nrcalong_rate.fits",
    ]
    for rate_file in rate_files:
        rtdata.get_data(rate_file)
        args = ["calwebb_image2", rtdata.input]
        Step.from_cmdline(args)

    # Get the level3 association json file (though not its members) and run
    # image3 pipeline on all _cal files listed in association
    rtdata.get_data("nircam/image/jw01538-o046_20230331t102920_image3_00009_asn.json")
    args = [
        "calwebb_image3",
        rtdata.input,
        "--steps.tweakreg.save_results=True",
        "--steps.source_catalog.snr_threshold=20",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)


def test_log_tracked_resources_detector1(log_tracked_resources, run_detector1pipeline_with_sirs):
    log_tracked_resources()


def test_log_tracked_resources_image2(log_tracked_resources, run_image2pipeline):
    log_tracked_resources()


def test_log_tracked_resources_image3(log_tracked_resources, run_image3pipeline):
    log_tracked_resources()


def test_nircam_image_sirs(run_detector1pipeline_with_sirs, rtdata_module, fitsdiff_default_kwargs):
    """Regression test of detector1 and image2 pipelines performed on NIRCam data."""
    rtdata = rtdata_module
    rtdata.input = "jw01345001001_10201_00001_nrca3_uncal.fits"
    output = "jw01345001001_10201_00001_nrca3_sirs_refpix.fits"
    rtdata.output = output
    rtdata.get_truth(
        "truth/test_nircam_image_stages/jw01345001001_10201_00001_nrca3_sirs_refpix.fits"
    )

    fitsdiff_default_kwargs["rtol"] = 5e-5
    fitsdiff_default_kwargs["atol"] = 1e-4

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "suffix",
    [
        "dq_init",
        "saturation",
        "superbias",
        "refpix",
        "linearity",
        "dark_current",
        "jump",
        "rate",
        "flat_field",
        "cal",
        "i2d",
    ],
)
def test_nircam_image_stages12(run_image2pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 and image2 pipelines performed on NIRCam data."""
    rtdata = rtdata_module
    rtdata.input = "jw01538046001_03105_00001_nrcalong_uncal.fits"
    output = "jw01538046001_03105_00001_nrcalong_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "i2d":
        fitsdiff_default_kwargs["rtol"] = 5e-5
        fitsdiff_default_kwargs["atol"] = 1e-4

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nircam_image_stage2_wcs(run_image2pipeline, rtdata_module):
    """Test that WCS object works as expected"""
    rtdata = rtdata_module
    rtdata.input = "jw01538046001_03105_00001_nrcalong_uncal.fits"
    output = "jw01538046001_03105_00001_nrcalong_assign_wcs.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    with datamodels.open(rtdata.output) as model, datamodels.open(rtdata.truth) as model_truth:
        grid = grid_from_bounding_box(model.meta.wcs.bounding_box)

        ra, dec = model.meta.wcs(*grid)
        ra_truth, dec_truth = model_truth.meta.wcs(*grid)

        assert_allclose(ra, ra_truth)
        assert_allclose(dec, dec_truth)


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
            if "jw01538046001_0310f_00001" not in model.meta.filename:
                assert "v2v3corr" in model.meta.wcs.available_frames


@pytest.mark.parametrize("suffix", ["i2d"])
def test_nircam_image_stage3(run_image3pipeline, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test that resampled i2d looks good for NIRCam imaging"""
    rtdata = rtdata_module
    rtdata.input = "jw01538-o046_20230331t102920_image3_00009_asn.json"
    output = f"jw01538-o046_t024_nircam_clear-f444w_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 2e-4

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nircam_image_stage3_catalog(run_image3pipeline, rtdata_module, diff_astropy_tables):
    rtdata = rtdata_module
    rtdata.input = "jw01538-o046_20230331t102920_image3_00009_asn.json"
    output = "jw01538-o046_t024_nircam_clear-f444w_cat.ecsv"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    assert diff_astropy_tables(rtdata.output, rtdata.truth, rtol=1e-4, atol=1e-5)


def test_nircam_image_stage3_segm(run_image3pipeline, rtdata_module, fitsdiff_default_kwargs):
    """Test that segmentation map looks good for NIRCam imaging"""
    rtdata = rtdata_module
    rtdata.input = "jw01538-o046_20230331t102920_image3_00009_asn.json"
    output = "jw01538-o046_t024_nircam_clear-f444w_segm.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_nircam_frame_averaged_darks(rtdata, fitsdiff_default_kwargs):
    """Test optional frame-averaged darks output from DarkCurrentStep"""
    rtdata.get_data("nircam/image/jw01205015001_03101_00001_nrcb1_ramp.fits")

    dark_file = "jw01205015001_03101_00001_nrcb1_frame_averaged_dark.fits"
    args = [
        "jwst.dark_current.DarkCurrentStep",
        rtdata.input,
        f"--dark_output={dark_file}",
    ]
    Step.from_cmdline(args)
    rtdata.output = dark_file

    rtdata.get_truth(f"truth/test_nircam_image/{dark_file}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_imaging_distortion(rtdata, fitsdiff_default_kwargs):
    """Verify that the distortion correction round trips."""
    rtdata.get_data("nircam/image/jw01538046002_02103_00002_nrca1_cal.fits")
    model = datamodels.open("jw01538046002_02103_00002_nrca1_cal.fits")
    wcsobj = model.meta.wcs
    sky_to_detector = wcsobj.get_transform("world", "detector")
    detector_to_sky = wcsobj.get_transform("detector", "world")

    # we'll use the crpix as the simplest reference point
    ra = model.meta.wcsinfo.crval1
    dec = model.meta.wcsinfo.crval2

    x, y = sky_to_detector(ra, dec)
    raout, decout = detector_to_sky(x, y)

    assert_allclose(x, model.meta.wcsinfo.crpix1 - 1)
    assert_allclose(y, model.meta.wcsinfo.crpix2 - 1)
    assert_allclose(raout, ra)
    assert_allclose(decout, dec)


@pytest.mark.parametrize(
    "suffix",
    [
        "cfn_clean_flicker_noise",
        "mask",
        "flicker_bkg",
        "flicker_noise",
        "cfn_ramp",
        "cfn_rate",
        "cfn_rateints",
    ],
)
def test_nircam_image_detector1_with_clean_flicker_noise(
    run_detector1_with_clean_flicker_noise, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Test detector1 pipeline for NIRCam imaging data with noise cleaning."""
    rtdata = rtdata_module
    rtdata.input = "jw01538046001_03105_00001_nrcalong_uncal.fits"
    output = f"jw01538046001_03105_00001_nrcalong_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nircam_image_clean_flicker_noise/{output}")

    # Set tolerances so the file comparisons work across architectures
    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 1e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("suffix", ["likely_rate", "likely_rateints"])
def test_nircam_image_detector1_with_likelihood(
    run_detector1_with_likelihood, rtdata_module, fitsdiff_default_kwargs, suffix
):
    """Test detector1 pipeline for NIRCam imaging data with likelihood fitting."""
    rtdata = rtdata_module
    rtdata.input = "jw01345001001_10201_00001_nrca3_uncal.fits"
    output = f"jw01345001001_10201_00001_nrca3_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nircam_image_likelihood/{output}")

    # Set tolerances so the file comparisons work across architectures
    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 1e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
