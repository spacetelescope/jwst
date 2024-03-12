"""
Unit test for Residual Fringe Correction fitting of the background
"""
import pytest

from pathlib import Path

from jwst.residual_fringe import utils
from jwst.residual_fringe.utils import fit_residual_fringes_1d as rf1d
from numpy.testing import assert_allclose
from astropy.io import fits


def read_fit_column(file):
    """ This is really  a small regression test, testing that the background fitting is working """

    # Data was pulled out of an exposure by modifying residual_fringe.py to write out a column of data
    # The function we are testing is fit_1d_background_complex.

    file_dir = Path(__file__).parent.resolve()
    file_path = str(file_dir / file)

    with fits.open(file_path) as hdu:
        col_data = hdu[1].data
        col_weight = hdu[2].data
        col_wnum = hdu[3].data
        bg_fit = hdu[4].data
        store_freq = hdu[0].header['FFREQ']

        return col_data, col_weight, col_wnum, bg_fit, store_freq


@pytest.mark.parametrize("file", ['good_col.fits', 'edge_col.fits'])
def test_background_fit(file):
    """ test fit_1d_background_complex"""

    (col_data, col_weight, col_wnum, bg_fit, store_freq) = read_fit_column(file)

    bg_fit2, _ = utils.fit_1d_background_complex(col_data, col_weight,
                                                 col_wnum, ffreq=store_freq)

    assert_allclose(bg_fit, bg_fit2, atol=0.001)

def read_spectrum(file):
    """ Read data from a small spectrum file for test_rf1d  """

    file_dir = Path(__file__).parent.resolve()
    file_path = str(file_dir / file)

    with fits.open(file_path) as hdu:
        wave = hdu[1].data
        influx = hdu[2].data
        outflux = hdu[3].data

        return wave, influx, outflux

@pytest.mark.parametrize("file", ['testspec.fits'])
def test_rf1d(file):
    """ Test the performance of the 1d residual defringe routine
        Do this using a Ch2C spectrum taken from observations of bright star 16 CygB
        We will test calling the 1d defringer directly.
    """

    (wave, influx, outflux) = read_spectrum(file)
    outflux2 = rf1d(influx, wave, channel=2)

    assert_allclose(outflux, outflux2, atol=0.001)
