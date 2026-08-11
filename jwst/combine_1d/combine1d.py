import logging

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer

log = logging.getLogger(__name__)

__all__ = [
    "count_input",
    "compute_output_wl",
    "check_exptime",
    "combine_1d_spectra",
    "check_monotonic",
]


def count_input(input_spectra):
    """
    Determine the number of input spectra that cover each wavelength.

    For any given input spectrum, the array of wavelengths gives the
    wavelengths at the centers of the pixels.  In this context, the
    expression that an input spectrum "covers" some particular wavelength
    (typically not one of the elements in the input spectrum's array of
    wavelengths) means that this wavelength is within the interval between
    the left edge of the first pixel and the right edge of the last pixel
    of the input spectrum.

    Parameters
    ----------
    input_spectra : list
        List of input spectra.

    Returns
    -------
    wl : ndarray
        Sorted list of all the wavelengths in all the input spectra.

    n_input_spectra : ndarray, 1-D
        For each element of ``wl``, the corresponding element of
        ``n_input_spectra`` is the number of input spectra that cover the
        wavelength in ``wl``.
    """
    # Create an array with all the input wavelengths (i.e. the union
    # of the input wavelengths).
    wl = None
    for in_spec in input_spectra:
        input_wl = in_spec.wavelength
        # only include spectra that have more than 1 data point
        if len(input_wl) > 1:
            if wl is None:
                wl = input_wl.copy()
            else:
                wl = np.hstack((input_wl, wl))
    wl.sort()
    nwl = len(wl)

    # n_input_spectra will be the number of input spectra that cover the
    # corresponding wavelength in wl.
    n_input_spectra = np.zeros(nwl, dtype=np.int64)
    for in_spec in input_spectra:
        input_wl = in_spec.wavelength

        # Check for degenerate spectrum. Skip with log.
        if len(input_wl) < 2:
            log.warning(f"Spectrum {in_spec} is degenerate with length {len(input_wl)}")
            log.warning("Skipping...")
            continue

        # wl0 and wl1 will be about a half pixel wider on either side
        # of the wavelength range for the current input spectrum.
        if input_wl[1] > input_wl[0]:  # wavelengths are increasing
            wl0 = input_wl[0] - 0.5 * (input_wl[1] - input_wl[0])
            wl1 = input_wl[-1] + 0.5 * (input_wl[-1] - input_wl[-2])
        elif input_wl[1] < input_wl[0]:  # wavelengths are decreasing
            wl0 = input_wl[-1] - 0.5 * (input_wl[-2] - input_wl[-1])
            wl1 = input_wl[0] + 0.5 * (input_wl[0] - input_wl[1])
        else:
            log.warning(f"Spectrum {in_spec} has a monotonic wavelength solution.")
            log.warning("Skipping...")
            continue
        temp = np.where(wl >= wl0, 1, 0)
        temp = np.where(wl >= wl1, 0, temp)
        n_input_spectra += temp
        del temp

    # This shouldn't happen.
    if np.any(n_input_spectra <= 0.0):
        raise RuntimeError("Problem with input wavelengths.")

    return wl, n_input_spectra


def compute_output_wl(wl, n_input_spectra):
    """
    Compute output wavelengths.

    In summary, the output wavelengths are computed by binning the
    input wavelengths in groups of the number of overlapping spectra.
    The merged and sorted input wavelengths are in ``wl``.

    If the wavelength arrays of the input spectra are nearly all the
    same, the values in wl will be in clumps of N nearly identical
    values (for N input files).  We want to average the wavelengths
    in those clumps to get the output wavelengths.  However, if the
    wavelength scales of the input files differ slightly, the values
    in ``wl`` may clump in some regions but not in others.  The algorithm
    below tries to find clumps; a clump is identified by having a
    small standard deviation of all the wavelengths in a slice of ``wl``,
    compared with neighboring slices.  It is intended that this should
    not find clumps where they're not present, i.e., in regions where
    the wavelengths of the input spectra are not well aligned with each
    other.  These regions will be gaps in the initial determination
    of output wavelengths based on clumps.  In such regions, output
    wavelengths will then be computed by just averaging groups of
    input wavelengths to fill in the gaps.

    Parameters
    ----------
    wl : 1-D array
        An array containing all the wavelengths from all the input
        spectra, sorted in increasing order.

    n_input_spectra : 1-D array
        An integer array of the same length as ``wl``.  For a given
        array index ``k`` (for example), ``n_input_spectra[k]`` is the number of
        input spectra that cover wavelength ``wl[k]``.

    Returns
    -------
    wavelength :  1-D array
        Array of wavelengths for the output spectrum.
    """
    nwl = len(wl)

    # sigma is an array of the standard deviation at each element
    # of wl, over n_input_spectra elements.  A small value implies that
    # there's a clump, i.e. several elements of wl with nearly the
    # same wavelength.
    sigma = np.zeros(nwl, dtype=np.float64) + 9999.0

    # mean_wl is the mean wavelength over the same slice of wl that we
    # used to compute sigma.  If sigma is small enough that it looks
    # as if there's a clump, we'll copy the mean_wl value to temp_wl
    # to be one element of the output wavelengths.
    mean_wl = np.zeros(nwl, dtype=np.float64) - 99.0

    # temp_wl has the same number of elements as wl, but we expect the
    # array of output wavelengths to be significantly smaller, so
    # temp_wl is initialized to a negative value as a flag.  Positive
    # elements will be copied to the array of output wavelengths.
    temp_wl = np.zeros(nwl, dtype=np.float64) - 99.0

    for k in range(nwl):
        n = n_input_spectra[k]
        if n == 1:
            sigma[k] = 0.0
            mean_wl[k] = wl[k]
            temp_wl[k] = mean_wl[k]
        else:
            k1 = k + n // 2 + 1
            k0 = k1 - n
            if k0 >= 0 and k1 <= nwl:
                sigma[k] = wl[k0:k1].std()
                mean_wl[k] = wl[k0:k1].mean()
                if sigma[k] == 0.0:
                    temp_wl[k] = mean_wl[k]

    cutoff = 0.8
    for k in range(nwl):
        # If sigma[k] equals 0, temp_wl has already been assigned.
        if sigma[k] > 0.0:
            if k == 0:
                if sigma[k] < cutoff * sigma[1]:
                    temp_wl[k] = mean_wl[k]
            elif k == nwl - 1:
                if sigma[k] < cutoff * sigma[nwl - 2]:
                    temp_wl[k] = mean_wl[k]
            else:
                if sigma[k] < cutoff * (sigma[k - 1] + sigma[k + 1]) / 2.0:
                    temp_wl[k] = mean_wl[k]

    # Fill gaps in the output wavelengths by taking averages of the
    # input wavelengths.  If there are n overlapping input spectra,
    # average a block of n elements of `wl`.
    done = False
    i = 0
    while not done:
        if i >= nwl:
            done = True
        else:
            n = n_input_spectra[i]
            middle = n // 2
            if i + middle < nwl:
                n = n_input_spectra[i + middle]
                if i + n < nwl:
                    # The nominal range is i - 1 to i + n (inclusive) for
                    # checking whether an output wavelength has already
                    # been assigned.
                    low = max(i - 1, 0)
                    high = min(i + n + 1, nwl - 1)
                    if temp_wl[low:high].max() <= 0.0:
                        temp_wl[i] = wl[i : i + n].mean()
            i += n

    return temp_wl[np.where(temp_wl > 0.0)].copy()


def check_exptime(exptime_key):
    """
    Check exptime_key for validity.

    This function checks ``exptime_key``.  If it is valid, the corresponding
    value used by the metadata interface will be returned.  This will be
    either "integration_time" or "exposure_time".  If it is invalid,
    "unit weight" will be returned (meaning that a weight of 1 will be
    used when averaging spectra), and a warning will be logged.

    Parameters
    ----------
    exptime_key : str
        A keyword or string indicating what value (integration time or
        exposure time) should be used as a weight when combing spectra.

    Returns
    -------
    exptime_key : str
        The value will be either "integration_time", "exposure_time",
        or "unit_weight".
    """
    exptime_lwr = exptime_key.lower()
    if exptime_lwr.startswith("integration") or exptime_lwr == "effinttm":
        exptime_key = "integration_time"
        log.info("Using integration time as the weight.")
    elif exptime_lwr.startswith("exposure") or exptime_lwr == "effexptm":
        exptime_key = "exposure_time"
        log.info("Using exposure time as the weight.")
    elif exptime_lwr == "unit_weight" or exptime_lwr == "unit weight":
        exptime_key = "unit_weight"
        log.info("Using weight = 1.")
    else:
        log.warning(f"Don't understand exptime_key = '{exptime_key}'; using unit weight.")
        log.info("The options for exptime_key are:")
        log.info("  integration_time, effinttm, exposure_time, effexptm, unit_weight, unit weight")
        exptime_key = "unit_weight"

    return exptime_key


def combine_1d_spectra(input_model, exptime_key, sigma_clip=None):
    """
    Combine the input spectra.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input spectra.  This will likely be a
        `~jwst.datamodels.container.ModelContainer` object,
        but may also be a multi-spectrum model, such as
        `~stdatamodels.jwst.datamodels.MultiSpecModel` or
        `~stdatamodels.jwst.datamodels.TSOMultiSpecModel`.
        Input spectra may have different spectral orders
        or wavelengths but should all share the same target.
        May be updated in place if processing is skipped.
    exptime_key : str
        A string identifying which keyword to use to get the exposure time,
        which is used as a weight when combining spectra.  The value should
        be one of:  "exposure_time" (the default), "integration_time",
        or "unit_weight".

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.MultiCombinedSpecModel`
        A combined spectra datamodel.
    """
    log.debug(f"Using exptime_key = {exptime_key}")

    exptime_key = check_exptime(exptime_key)

    input_spectra = {}
    output_spectra = {}
    if isinstance(input_model, ModelContainer):
        for ms in input_model:
            _read_input_spectra(ms, exptime_key, input_spectra)
    else:
        _read_input_spectra(input_model, exptime_key, input_spectra)

    if len(input_spectra) == 0:
        log.error("No valid input spectra found for source. Skipping.")
        input_model.meta.cal_step.combine_1d = "SKIPPED"
        return input_model

    for order in input_spectra:
        output_spectra[order] = OutputSpectrumModel()
        output_spectra[order].assign_wavelengths(input_spectra[order])
        output_spectra[order].accumulate_sums(input_spectra[order], sigma_clip=sigma_clip)

    output_model = datamodels.MultiCombinedSpecModel()

    for order in output_spectra:
        output_order = output_spectra[order].create_output_data()
        output_order.spectral_order = order
        for attr in SPECMETA_ATTRIBUTES:
            setattr(output_order, attr, getattr(input_spectra[order][0], attr))
        output_model.spec.append(output_order)

    # Copy one of the input headers to output.
    if isinstance(input_model, ModelContainer):
        output_model.update(input_model[0], only="PRIMARY")
    else:
        output_model.update(input_model, only="PRIMARY")

    # Looks clunky, but need an output_spec instance to copy wcs
    output_model.meta.wcs = output_spectra[list(output_spectra)[0]].wcs
    output_model.meta.cal_step.combine_1d = "COMPLETE"

    for order in input_spectra:
        for in_spec in input_spectra[order]:
            in_spec.close()

    for order in output_spectra:
        output_spectra[order].close()

    return output_model


def check_monotonic(wavelength):
    """
    Check if wavelength array is strictly monotonic (purely increasing or decreasing).

    Parameters
    ----------
    wavelength : list or array
        An array of wavelengths

    Returns
    -------
    bool
        True if the array is strictly increasing or strictly decreasing,
        False otherwise. Note that duplicates will return False.
    """
    if len(wavelength) < 2:
        return True

    is_increasing = all(np.diff(wavelength) > 0)
    is_decreasing = all(np.diff(wavelength) < 0)

    return is_increasing or is_decreasing
