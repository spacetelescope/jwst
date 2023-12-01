import numpy as np
import logging
from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(sp_leak_ref, ch1b, ch3a):
    """
    Short Summary
    -------------
    Using the spectral leak reference correction and spectrum containing CH1 B
    correct the CH 3 A spectrum.

    Parameters
    ----------
    sp_leak_ref : string
        name of the spectral leak reference file defined by datamodels.MirMrsPtCorrModel

    ch1b : numpy array
        Input channel 1 B spectrum
    ch3a : numpy array
        Input channel 3 A spectrum

    Returns
    -------
    output_model : ~jwst.datamodels.RampModel
        Spectral leak corrected science data
    """

    wave1b = ch1b.spec[0].spec_table.WAVELENGTH
    spec1b = ch1b.spec[0].spec_table.FLUX

    wave3a = ch3a.spec[0].spec_table.WAVELENGTH
    spec3a = ch3a.spec[0].spec_table.FLUX

    leak_ref = datamodels.MirMrsPtCorrModel(sp_leak_ref)
    leak_wave = leak_ref.leakcor_table.wavelength
    leak_percent = leak_ref.leakcor_table.frac_leak

    # Spectral leak vector

    leak = np.interp(wave3a, 2*wave1b, spec1b) * np.interp(wave3a, leak_wave, leak_percent)
    # Factor of 2 in 2*wavelb is  necessary because we're dealing with second order light.  Without this the interpolation simply
    # uses the last point in the spectrum, which is close enough to correct as won't be noticeable for BAND spectra,
    # but for CH spectra can be quite different.

    # Correct the data
    spec3a_corr = spec3a - leak
    return spec3a_corr
