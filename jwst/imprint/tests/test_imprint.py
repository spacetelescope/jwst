"""
Unit tests for imprint correction
"""

from jwst.datamodels import ImageModel
from jwst.imprint import ImprintStep
import numpy as np
import pytest


def test_step(make_imagemodel):
    """Assert that the results should be all zeros.
    """
    im = make_imagemodel(10, 10)
    pos_no = []
    imprint = []
    pos_no.append(0)
    imprint.append(im)
    result = ImprintStep.call(im, imprint, pos_no)

    assert result.meta.cal_step.imprint == 'COMPLETE'
    assert result.data.sum() == 0


@pytest.fixture(scope='function')
def make_imagemodel():
    '''Image model for testing'''
    def _im(ysize, xsize):
        # create the data arrays
        im = ImageModel((ysize, xsize))
        im.data = np.random.rand(ysize, xsize)

        return im

    return _im
