
import logging
import numpy as np

from . import utils
from . import nrm_consts

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def get_webbpsf_filter(filter_model, specbin=None, trim=False):
    """
    Short Summary
    -------------
    Load the filter bandpass data using files from WebbPSF

    Parameters
    ----------
    filter_model: filter model object
        filter throughput reference data

    specbin: integer
        factor to bin spectrum down by

    trim: boolean
        trims below 0.8 and above 1.2 lambda_c (dg - is this right?)

    Returns
    -------
    spec - 2D float array
        filter bandpass data: [weight, wavelength_in_microns]
    """
    W = 1  # remove confusion - wavelength index
    T = 0  # remove confusion - trans index after reading in...

    nwave = len(filter_model.filter_table['wavelength'])
    tmp_array = np.zeros((nwave, 2))

    # input in angst _ANGSTROM = 1.0e-10
    tmp_array[:, W] = filter_model.filter_table['wavelength'] * 1.0e-10

    # weights (peak unity, unnormalized sum)
    tmp_array[:, T] = filter_model.filter_table['throughput']

    # remove leading and trailing throughput lines with 'flag' array of indices
    flag = np.where(tmp_array[:, T] != 0)[0]
    spec = tmp_array[flag, :]

    # rebin as desired - fewer wavelengths for debugging quickly
    if specbin is not None:
        smallshape = spec.shape[0] // specbin

        log.debug(' bin filter data by %d from %d to %d', specbin,
                  spec.shape[0], smallshape)

        spec = spec[:smallshape * specbin, :]  # clip trailing
        spec = utils.krebin(spec, (smallshape, 2))
        spec[:, W] = spec[:, W] / float(specbin)  # krebin added up waves
        spec[:, T] = spec[:, T] / float(specbin)  # krebin added up trans too
        log.debug(' post-bin shape %s', spec.shape)

    if trim:
        log.debug(' TRIMing filter data ...')
        wl = spec[:, W].copy()
        tr = spec[:, T].copy()
        idx = np.where((wl > (1.0 - 0.5 * trim[1]) * trim[0])
                       & (wl < (1.0 + 0.5 * trim[1]) * trim[0]))
        wl = wl[idx]
        tr = tr[idx]
        spec = np.zeros((len(idx[0]), 2))
        spec[:, 1] = wl
        spec[:, 0] = tr
        log.debug(' post-trim shape %s', spec.shape)

    log.debug(' final filter shape %s', spec.shape)
    log.debug(' %d filter samples between %.3f and %.3f um',
              len(spec[:, 0]), spec[0, W] / nrm_consts.um_,
              spec[-1, W] / nrm_consts.um_)

    return spec
