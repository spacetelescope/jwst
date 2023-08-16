import numpy as np


WFSS_EXPTYPES = ['NIS_WFSS', 'NRC_WFSS', 'NRC_GRISM', 'NRC_TSGRISM']


def get_wavelengths(model, exp_type="", order=None, use_wavecorr=None):
    """Read or compute wavelengths.

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel`
        The input science data, or a slit from a
        `~jwst.datamodels.MultiSlitModel`.

    exp_type : str
        The exposure type.  This is only needed to check whether the input
        data are WFSS.

    order : int
        Spectral order number, for NIRISS SOSS only.

    use_wavecorr : bool
        Use the corrected wavelengths in the wavelength attribute or
        recompute uncorrected wavelengths from the WCS?

    Returns
    -------
    wl_array : 2-D ndarray
        An array of wavelengths corresponding to the data in `model`.
    """

    # Use the existing wavelength array, if there is one
    if hasattr(model, "wavelength"):
        wl_array = model.wavelength.copy()
        got_wavelength = True                   # may be reset below
    else:
        wl_array = None
    if (wl_array is None or len(wl_array) == 0 or np.nanmin(wl_array) == 0.
            and np.nanmax(wl_array) == 0.):
        got_wavelength = False
        wl_array = None

    # If we've been asked to use the corrected wavelengths stored in
    # the wavelength array, return those wavelengths. Otherwise, the
    # results computed from the WCS object (below) will be returned.
    if use_wavecorr is not None:
        if use_wavecorr:
            return wl_array
        else:
            got_wavelength = False  # force wl computation below

    # Evaluate the WCS on the grid of pixel indexes, capturing only the
    # resulting wavelength values
    shape = model.data.shape
    grid = np.indices(shape[-2:], dtype=np.float64)

    # If no existing wavelength array, compute one
    if hasattr(model.meta, "wcs") and not got_wavelength:
        # Set up an appropriate WCS object
        if hasattr(model.meta, "exposure") and model.meta.exposure.type == "NIS_SOSS":
            wl_array = model.meta.wcs(grid[1], grid[0], order)[2]
            return wl_array

        wcs = model.meta.wcs

        if exp_type in WFSS_EXPTYPES:
            # We currently have to loop over pixels for WFSS data.
            wl_array = np.zeros(shape[-2:], dtype=np.float64)
            for j in range(shape[-2]):
                for i in range(shape[-1]):
                    # Keep wavelength; ignore RA and Dec
                    wl_array[..., j, i] = wcs(i, j)[2]
        else:
            wl_array = wcs(grid[1], grid[0])[2]

    return wl_array
