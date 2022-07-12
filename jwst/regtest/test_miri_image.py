import pytest
from astropy.io.fits.diff import FITSDiff
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box
from jwst.stpipe import Step
from jwst import datamodels


@pytest.fixture(scope="module")
def run_detector1(rtdata_module):
    """Run detector1 pipeline on MIRI imaging data."""
    rtdata = rtdata_module
    rtdata.get_data("miri/image/jw01024001001_04101_00001_mirimage_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files
    args = ["jwst.pipeline.Detector1Pipeline", rtdata.input,
            "--save_calibrated_ramp=True",
            "--steps.dq_init.save_results=True",
            "--steps.saturation.save_results=True",
            "--steps.firstframe.save_results=True",
            "--steps.lastframe.save_results=True",
            "--steps.reset.save_results=True",
            "--steps.linearity.save_results=True",
            "--steps.rscd.save_results=True",
            "--steps.dark_current.save_results=True",
            "--steps.refpix.save_results=True",
            ]
    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_image2(run_detector1, rtdata_module):
    """Run image2 pipeline on the _rate file, saving intermediate products"""
    rtdata = rtdata_module
    rtdata.input = 'jw01024001001_04101_00001_mirimage_rate.fits'
    args = ["jwst.pipeline.Image2Pipeline", rtdata.input,
            "--steps.assign_wcs.save_results=True",
            "--steps.flat_field.save_results=True"
            ]
    Step.from_cmdline(args)

    # Grab rest of _rate files for the asn and run image2 pipeline on each to
    # produce fresh _cal files for the image3 pipeline.  We won't check these
    # or look at intermediate products, and skip resample (don't need i2d image)
    rate_files = [
        "miri/image/jw01024001001_04101_00002_mirimage_rate.fits",
        "miri/image/jw01024001001_04101_00003_mirimage_rate.fits",
        "miri/image/jw01024001001_04101_00004_mirimage_rate.fits",
    ]
    for rate_file in rate_files:
        rtdata.get_data(rate_file)
        args = ["jwst.pipeline.Image2Pipeline", rtdata.input,
                "--steps.resample.skip=True"]
        Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_image3(run_image2, rtdata_module):
    """Get the level3 association json file (though not its members) and run
    image3 pipeline on all _cal files listed in association"""
    rtdata = rtdata_module
    rtdata.get_data("miri/image/jw01024-o001_20220501t155404_image3_001_asn.json")
    args = ["jwst.pipeline.Image3Pipeline", rtdata.input]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["dq_init", "saturation", "firstframe", "lastframe", "reset",
                                    "linearity", "rscd", "dark_current", "refpix", "ramp", "rate",
                                    "rateints"])
def test_miri_image_detector1(run_detector1, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 pipeline performed on MIRI imaging data."""
    _assert_is_same(rtdata_module, fitsdiff_default_kwargs, suffix)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["assign_wcs", "flat_field", "cal", "i2d"])
def test_miri_image_image2(run_image2, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of image2 pipeline performed on MIRI imaging data."""
    _assert_is_same(rtdata_module, fitsdiff_default_kwargs, suffix)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["o001_crf"])
def test_miri_image_image3(run_image3, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of image3 pipeline performed on MIRI imaging data."""
    _assert_is_same(rtdata_module, fitsdiff_default_kwargs, suffix)


def _assert_is_same(rtdata_module, fitsdiff_default_kwargs, suffix):
    """Assertion helper for the above tests"""
    rtdata = rtdata_module
    rtdata.input = "jw01024001001_04101_00001_mirimage_uncal.fits"
    output = f"jw01024001001_04101_00001_mirimage_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/test_miri_image_stages/{output}")

    # Set tolerances so the crf, rscd and rateints file comparisons work across
    # architectures
    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 1e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_image3_i2d(run_image3, rtdata_module, fitsdiff_default_kwargs):
    rtdata = rtdata_module
    rtdata.input = "jw01024-o001_20220501t155404_image3_001_asn.json"
    rtdata.output = "jw01024-o001_t002_miri_f770w_i2d.fits"
    rtdata.get_truth("truth/test_miri_image_stages/jw01024-o001_t002_miri_f770w_i2d.fits")

    fitsdiff_default_kwargs["rtol"] = 1e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_image3_catalog(run_image3, rtdata_module, diff_astropy_tables):
    rtdata = rtdata_module
    rtdata.input = "jw01024-o001_20220501t155404_image3_001_asn.json"
    rtdata.output = "jw01024-o001_t002_miri_f770w_cat.ecsv"
    rtdata.get_truth("truth/test_miri_image_stages/jw01024-o001_t002_miri_f770w_cat.ecsv")

    assert diff_astropy_tables(rtdata.output, rtdata.truth, rtol=1e-4)


@pytest.mark.bigdata
def test_miri_image_wcs(run_image2, rtdata_module, fitsdiff_default_kwargs):
    rtdata = rtdata_module

    # get input assign_wcs and truth file
    output = "jw01024001001_04101_00001_mirimage_assign_wcs.fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_miri_image_stages/" + output)

    # Open the output and truth file
    with datamodels.open(rtdata.output) as im, datamodels.open(rtdata.truth) as im_truth:
        x, y = grid_from_bounding_box(im.meta.wcs.bounding_box)
        ra, dec = im.meta.wcs(x, y)
        ratruth, dectruth = im_truth.meta.wcs(x, y)
        assert_allclose(ra, ratruth)
        assert_allclose(dec, dectruth)

        # Test the inverse transform
        xtest, ytest = im.meta.wcs.backward_transform(ra, dec)
        xtruth, ytruth = im_truth.meta.wcs.backward_transform(ratruth, dectruth)
        assert_allclose(xtest, xtruth)
        assert_allclose(ytest, ytruth)
