import numpy as np
from numpy.testing import assert_allclose

from ... master_background import create_master_bkg
from ...combine_1d import combine1d


def test_to_pixel():
    wave = np.linspace(.3, 6.8, 10)
    flux = np.arange(10, dtype=np.float64) + 100
    bkg = create_master_bkg.create_background(wave, flux)
    inp = combine1d.InputSpectrumModel(bkg, bkg.spec[0],
                                       exptime_key='unit_weight',
                                       background=False)
    outm = combine1d.OutputSpectrumModel()
    outm.assign_wavelengths([inp])
    wave = [.3, 6.8, 2.2, .1]
    res_to_pixel = outm.to_pixel(wave)
    res_wcs = outm.wcs.invert(wave, wave, wave)
    assert_allclose(res_to_pixel, res_wcs)
    

def test_combine_1d():
    wave1 = np.linspace(1.2, 5.6, 10)
    wave2 = np.linspace(.3, 6.8, 10)
    flux = np.arange(10, dtype=np.float64) + 100

    bkg1 = create_master_bkg.create_background(wave1, flux)
    bkg2 = create_master_bkg.create_background(wave2, flux)

    inp2 = combine1d.InputSpectrumModel(bkg2, bkg2.spec[0],
                                        exptime_key='unit_weight',
                                        background=False)
    inp1 = combine1d.InputSpectrumModel(bkg1, bkg1.spec[0],
                                        exptime_key='unit_weight',
                                        background=False)
    outm = combine1d.OutputSpectrumModel()
    outm.assign_wavelengths([inp1, inp2])
    outm.accumulate_sums([inp1, inp2], 'linear')
    outm.compute_combination()
