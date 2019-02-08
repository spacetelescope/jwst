import logging

import numpy as np

from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def create_background(wavelength, flux):
    """Create a 1-D spectrum table as a MultiSpecModel.

    This is the syntax for accessing the data in the columns:
    wavelength = output_model.spec[0].spec_table['wavelength']
    background = output_model.spec[0].spec_table['flux']

    Parameters
    ----------
    wavelength : 1-D ndarray
        Array of wavelengths, in micrometers.

    flux : 1-D ndarray
        Array of background fluxes.

    Returns
    -------
    output_model : `~jwst.datamodels.MultiSlitModel`, or None
        A data model containing the 1-D background spectrum.  This can be
        written to disk by calling:

        output_model.save(<filename>)
    """

    wl_shape = wavelength.shape
    flux_shape = flux.shape
    bad = False
    if len(wl_shape) > 1:
        bad = True
        log.error("The wavelength array has shape {}; expected "
              "a 1-D array".format(wl_shape))
    if len(flux_shape) > 1:
        bad = True
        log.error("The background flux array has shape {}; expected "
              "a 1-D array".format(flux_shape))
    if bad:
        return None

    if wl_shape[0] != flux_shape[0]:
        log.error("wavelength array has length {}, "
                  "but background flux array has length {}."
                  .format(wl_shape[0], flux_shape[0]))
        log.error("The arrays must be the same size.")
        return None

    # Create arrays for columns that we won't need.
    dummy = np.zeros(wl_shape[0], dtype=np.float64)
    dq = np.zeros(wl_shape[0], dtype=np.int32)

    output_model = datamodels.MultiSpecModel()

    spec_dtype = datamodels.SpecModel().spec_table.dtype

    # xxx Include one more argument at the end, after column NPIXELS has
    # xxx been added to the x1d table.  The new argument should be float64
    # xxx and with values of 1 (i.e. not dummy).
    otab = np.array(list(zip(wavelength, flux,
                             dummy, dq, dummy, dummy, dummy, dummy)),
                    dtype=spec_dtype)

    spec = datamodels.SpecModel(spec_table=otab)
    output_model.spec.append(spec)

    return output_model
