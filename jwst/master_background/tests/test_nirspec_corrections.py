"""
Unit tests for master background NIRSpec corrections
"""
import numpy as np

from stdatamodels.jwst import datamodels

from jwst.master_background.nirspec_utils import (
    correct_nrs_ifu_bkg, correct_nrs_fs_bkg
)


def test_ifu_pathloss_existence():
    """Test the case where the input is missing a pathloss array"""

    input = datamodels.IFUImageModel((10, 10))
    result = correct_nrs_ifu_bkg(input)

    assert result == input


def test_ifu_correction():
    """Test application of IFU corrections"""

    data = np.ones((5, 5))
    pl_ps = 2.1 * data
    pl_un = data / 1.9
    input = datamodels.IFUImageModel(data=data,
                                     pathloss_point=pl_ps,
                                     pathloss_uniform=pl_un)

    corrected = input.data * pl_un / pl_ps
    result = correct_nrs_ifu_bkg(input)

    assert np.allclose(corrected, result.data, rtol=1.e-7)


def test_fs_correction():
    """Test application of FS corrections"""

    data = np.ones((5, 5))
    ff_ps = 1.5 * data
    ff_un = data / 1.2
    pl_ps = 2 * data
    pl_un = data / 2.1
    ph_ps = 1.1 * data
    ph_un = 1.23 * data
    input = datamodels.SlitModel(data=data,
                                 flatfield_point=ff_ps, flatfield_uniform=ff_un,
                                 pathloss_point=pl_ps, pathloss_uniform=pl_un,
                                 photom_point=ph_ps, photom_uniform=ph_un)

    corrected = input.data * (ff_un / ff_ps) * (pl_un / pl_ps) * (ph_ps / ph_un)
    result = correct_nrs_fs_bkg(input, primary_slit=True)

    assert np.allclose(corrected, result.data, rtol=1.e-7)
