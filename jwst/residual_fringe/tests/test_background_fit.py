"""Unit test for Residual Fringe Correction fitting of the background."""

import pytest
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from numpy.testing import assert_allclose

from jwst.residual_fringe import utils


def read_fit_column(file):
    """
    Read some small sample data for testing.

    Parameters
    ----------
    file : str
       File name, should be stored in the same directory
       as this file.

    Returns
    -------
    tuple
        Test data: col_data, col_weight, col_wnum, bg_fit, store_freq.
    """
    # Data was pulled out of an exposure by modifying
    # residual_fringe.py to write out a column of data
    # The function we are testing is fit_1d_background_complex.

    file_path = get_pkg_data_filename(f"data/{file}", package="jwst.residual_fringe.tests")

    with fits.open(file_path) as hdu:
        col_data = hdu[1].data
        col_weight = hdu[2].data
        col_wnum = hdu[3].data
        bg_fit = hdu[4].data
        store_freq = hdu[0].header["FFREQ"]

    return col_data, col_weight, col_wnum, bg_fit, store_freq


@pytest.mark.parametrize("file", ["good_col.fits", "edge_col.fits"])
def test_background_fit(file):
    (col_data, col_weight, col_wnum, bg_fit, store_freq) = read_fit_column(file)

    bg_fit2, _ = utils.fit_1d_background_complex(col_data, col_weight, col_wnum, ffreq=store_freq)

    assert_allclose(bg_fit, bg_fit2, atol=0.001)
