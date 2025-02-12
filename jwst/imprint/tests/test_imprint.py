"""Unit tests for imprint correction."""

from stdatamodels.jwst.datamodels import ImageModel

from jwst.imprint import ImprintStep
import numpy as np


def make_imagemodel(ysize, xsize, value=None):
    """Make an image model for testing."""
    im = ImageModel()
    if value is None:
        rng = np.random.default_rng(seed=42)
        im.data = rng.random((ysize, xsize))
    else:
        im.data = np.full((ysize, xsize), value)
    return im


def test_step():
    """Test that the results are all zeros when subtracting an image from itself."""
    im = make_imagemodel(10, 10)
    imprint = [im]
    result = ImprintStep.call(im, imprint)

    assert result.meta.cal_step.imprint == "COMPLETE"
    assert result.data.sum() == 0


def test_step_single_imprint():
    # Make science and imprint with mismatched background flags
    science = make_imagemodel(10, 10)
    science.meta.observation.bkgdtarg = False
    imprint = make_imagemodel(10, 10)
    imprint.meta.observation.bkgdtarg = True

    result = ImprintStep.call(science, [imprint])

    # for a single imprint, it's used anyway
    assert result.meta.cal_step.imprint == "COMPLETE"
    assert result.data.sum() == 0


def test_step_match_dither():
    # Make imprints with a range of dither positions
    science = make_imagemodel(10, 10, value=3.0)
    science.meta.dither.position_number = 2

    imprints = []
    for i in range(4):
        imprint = make_imagemodel(10, 10, value=i)
        imprint.meta.dither.position_number = i + 1
        imprints.append(imprint)

    result = ImprintStep.call(science, imprints)

    # The matching imprint is used (i=1, value=1.0)
    assert result.meta.cal_step.imprint == "COMPLETE"
    assert np.all(result.data == 2.0)


def test_step_match_background():
    # Make science and imprint with mismatched background flags
    science = make_imagemodel(10, 10, value=3.0)
    science.meta.observation.bkgdtarg = False
    imprint_bg = make_imagemodel(10, 10, value=2.0)
    imprint_bg.meta.observation.bkgdtarg = True
    imprint_sci = make_imagemodel(10, 10, value=1.0)
    imprint_sci.meta.observation.bkgdtarg = False

    result = ImprintStep.call(science, [imprint_bg, imprint_sci])

    # The matching imprint is used
    assert result.meta.cal_step.imprint == "COMPLETE"
    assert np.all(result.data == 2.0)


def test_step_match_background_mismatched_dither():
    # Make science and imprint with mismatched background flags
    # and dither positions
    science = make_imagemodel(10, 10, value=3.0)
    science.meta.observation.bkgdtarg = False
    science.meta.dither.position_number = 2
    imprint_bg = make_imagemodel(10, 10, value=2.0)
    imprint_bg.meta.observation.bkgdtarg = True
    imprint_bg.meta.dither.position_number = 1
    imprint_sci = make_imagemodel(10, 10, value=1.0)
    imprint_sci.meta.observation.bkgdtarg = False
    imprint_sci.meta.dither.position_number = 1

    result = ImprintStep.call(science, [imprint_bg, imprint_sci])

    # The matching imprint is used
    assert result.meta.cal_step.imprint == "COMPLETE"
    assert np.all(result.data == 2.0)


def test_step_no_match():
    science = make_imagemodel(10, 10, value=3.0)
    science.meta.dither.position_number = 5

    # Make imprints with a range of dither positions,
    # none of which match the science
    imprints = []
    for i in range(4):
        imprint = make_imagemodel(10, 10, value=i)
        imprint.meta.dither.position_number = i + 1
        imprints.append(imprint)

    result = ImprintStep.call(science, imprints)

    # No match is found
    assert result.meta.cal_step.imprint == "SKIPPED"
    assert np.all(result.data == 3.0)
