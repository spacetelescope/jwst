import numpy as np
import logging
import math
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
    sp_leak_ref: ~datamodels.MirMrsPtCorrModel
        spectral leak reference information

    ch1b: numpy array 

    ch3a: numpy array

    Returns
    -------
    output_model: ~jwst.datamodels.RampModel
        RSCD-corrected science data
    """

    wave1b = ch1b.spec[0].spec_table.WAVELENGTH
    spec1b = ch1b.spec[0].spec_table.FLUX

    wave3a = ch3a.spec[0].spec_table.WAVELENGTH
    spec3a = ch3a.spec[0].spec_table.FLUX

    leak_ref = datamodels.MirMrsPtCorrModel(sp_leak_ref)
    leak_wave = leak_ref.leakcor_table.wavelength
    leak_percent = leak_ref.leakcor_table.frac_leak

    # Spectral leak vector
    leak = np.interp(wave3a, wave1b, spec1b) * np.interp(wave3a, leak_wave, leak_percent)

    # Correct the data
    spec3a_corr = spec3a - leak
    return spec3a_corr
