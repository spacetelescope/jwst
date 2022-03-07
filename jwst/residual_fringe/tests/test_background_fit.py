"""
Unit test for Residual Fringe Correction fitting of the background
"""

from pathlib import Path

from jwst.residual_fringe import utils
from numpy.testing import assert_allclose
from astropy.io import fits


def read_fit_column(file):
    """ This is really  a small regression test, testing that the background fitting is working """

    # Data was pulled out of an exposure by modifying residual_fringe.py to write out a column of data
    # The function we are testing is fit_1d_background_complex.

    file_dir = Path(__file__).parent.resolve()
    file_path = str(file_dir / file)

    hdu = fits.open(file_path)
    col_data = hdu[1].data
    col_weight = hdu[2].data
    col_wnum = hdu[3].data
    bg_fit = hdu[4].data
    store_freq = hdu[0].header['FFREQ']
    bg_fit2, _ = utils.fit_1d_background_complex(col_data, col_weight,
                                                 col_wnum, ffreq=store_freq)

    assert_allclose(bg_fit, bg_fit2, atol=0.001)


def test_background_fit():
    """ test fit_1d_background_complex"""

    # test a redefined good column
    file1 = 'residual_col_633.fits'
    read_fit_column(file1)

    # test a redefined edge column (sometimes tricky to fit)
    file2 = 'residual_col_433.fits'
    read_fit_column(file2)
