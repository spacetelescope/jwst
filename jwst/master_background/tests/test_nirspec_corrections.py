"""
Unit tests for master background NIRSpec corrections
"""
import numpy as np

from jwst import datamodels
from ..nirspec_utils import correct_nrs_ifu_bkg


def test_ifu_pathloss_existence():
    """Test the case where the input is missing a pathloss array"""

    input = datamodels.IFUImageModel((10, 10))
    result = correct_nrs_ifu_bkg(input)

    assert (result == input)


def test_ifu_correction():
    """Test application of IFU corrections"""

    data = np.ones((5, 5))
    pl_ps = 2 * data
    pl_un = data / 2
    input = datamodels.IFUImageModel(data=data,
                                     pathloss_point=pl_ps,
                                     pathloss_uniform=pl_un)

    corrected = input.data * pl_ps / pl_un
    result = correct_nrs_ifu_bkg(input)

    assert np.allclose(corrected, result.data, rtol=1.e-10)
