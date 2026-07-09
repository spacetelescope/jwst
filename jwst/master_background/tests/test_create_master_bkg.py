"""Test master_background.create_master_bkg."""

import numpy as np

from jwst.master_background import create_master_bkg


def test_create_master_bkg(caplog):
    wl = np.arange(37, dtype=np.float64) * 0.1 + 3.0
    surf_bright = np.ones(37, dtype=np.float64) * 7

    # Bad: mismatch between array lengths
    model = create_master_bkg.create_background(wl, surf_bright[0:5])
    assert model is None
    assert "arrays must be the same size" in caplog.text
    caplog.clear()

    # Bad: wrong wavelength shape
    mock_wl = np.zeros((5, 7), dtype=np.float64) + 3
    model = create_master_bkg.create_background(mock_wl, surf_bright)
    assert model is None
    assert "wavelength array has shape (5, 7); expected a 1-D array" in caplog.text
    caplog.clear()

    # Bad: wrong surf_bright shape
    mock_sb = np.zeros((5, 7), dtype=np.float64) + 3
    model = create_master_bkg.create_background(wl, mock_sb)
    assert model is None
    assert "surf_bright array has shape (5, 7); expected a 1-D array" in caplog.text
    caplog.clear()

    # Good: all data match and have expected shapes
    model = create_master_bkg.create_background(wl, surf_bright)

    tab_wavelength = model.spec[0].spec_table["wavelength"]
    tab_surf_bright = model.spec[0].spec_table["surf_bright"]
    assert np.allclose(wl, tab_wavelength, rtol=1.0e-10)
    assert np.allclose(surf_bright, tab_surf_bright, rtol=1.0e-10)
