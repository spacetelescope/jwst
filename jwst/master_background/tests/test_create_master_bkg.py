"""
Test for master_background.create_master_bkg
"""
import numpy as np

from jwst.master_background import create_master_bkg

def test_create_master_bkg():

    wl = np.arange(37, dtype=np.float64) * 0.1 + 3.
    flux = np.ones(37, dtype=np.float64) * 7

    model = create_master_bkg.create_background(wl, flux[0:5])
    assert model is None

    dummy = np.zeros((5, 7), dtype=np.float64) + 3
    model = create_master_bkg.create_background(dummy, flux)
    assert model is None

    model = create_master_bkg.create_background(wl, flux)

    tab_wavelength = model.spec[0].spec_table['wavelength']
    tab_flux = model.spec[0].spec_table['flux']
    assert np.allclose(wl, tab_wavelength, rtol=1.e-10)
    assert np.allclose(flux, tab_flux, rtol=1.e-10)
