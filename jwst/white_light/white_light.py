"""Sum the flux over all wavelengths in each integration as a function of time for the target."""

import logging

import numpy as np
from collections import OrderedDict
from astropy.table import QTable
from astropy.time import Time, TimeDelta

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def white_light(input_model, min_wave=None, max_wave=None):
    """
    Compute the integrated flux over all wavelengths for a multi-integration extracted spectrum.

    Parameters
    ----------
    input_model : MultiSpecModel
        Datamodel containing the multi-integration data

    min_wave : float, optional
        Default wavelength minimum for integration.

    max_wave : float, optional
        Default wavelength maximum for integration.

    Returns
    -------
    tbl : astropy.table.table.QTable
        Table containing the integrated flux as a function of time.
    """
    ntables = len(input_model.spec)
    fluxsums = []

    # The input should contain one row per integration for each spectral
    # order.  NIRISS SOSS data can contain up to three orders.
    norders = 0  # number of different spectral orders
    sporders = []  # list of spectral order numbers
    ntables_order = []  # number of tables for each spectral order
    prev_spectral_order = -999
    for i in range(ntables):
        # The following assumes that all rows for a given spectral order
        # are contiguous.
        spectral_order = input_model.spec[i].spectral_order
        if spectral_order is None:
            norders = 1
            sporders = [0]
            ntables_order = [ntables]
            break
        if spectral_order != prev_spectral_order:
            sporders.append(spectral_order)
            prev_spectral_order = spectral_order
            ntables_order.append(1)
            norders = len(sporders)
        else:
            ntables_order[norders - 1] += 1

    log.debug(
        "norders = %d, sporders = %s, ntables_order = %s",
        norders,
        str(sporders),
        str(ntables_order),
    )

    # Create a wavelength mask, using cutoffs if specified, then
    # compute the flux sum for each integration in the input.
    low_cutoff = -1.0
    high_cutoff = 1.0e10
    if min_wave is not None:
        low_cutoff = min_wave
    if max_wave is not None:
        high_cutoff = max_wave

    for i in range(ntables):
        wave_mask = np.where(
            (input_model.spec[i].spec_table["WAVELENGTH"] >= low_cutoff)
            & (input_model.spec[i].spec_table["WAVELENGTH"] <= high_cutoff)
        )[0]

        fluxsums.append(np.nansum(input_model.spec[i].spec_table["FLUX"][wave_mask]))

    # Populate metadata for the output table
    tbl_meta = OrderedDict()
    tbl_meta["instrument"] = input_model.meta.instrument.name
    tbl_meta["detector"] = input_model.meta.instrument.detector
    tbl_meta["exp_type"] = input_model.meta.exposure.type
    tbl_meta["subarray"] = input_model.meta.subarray.name
    tbl_meta["filter"] = input_model.meta.instrument.filter
    tbl_meta["pupil"] = input_model.meta.instrument.pupil
    tbl_meta["target_name"] = input_model.meta.target.catalog_name

    # Create the output table
    tbl = QTable(meta=tbl_meta)

    if hasattr(input_model, "int_times") and input_model.int_times is not None:
        nrows = len(input_model.int_times)
    else:
        nrows = 0
    if nrows == 0:
        log.warning("There is no INT_TIMES table in the input file.")

    if nrows > 0:
        int_start = input_model.meta.exposure.integration_start  # one indexed
        if int_start is None:
            int_start = 1
            log.warning("INTSTART not found; assuming a value of %d", int_start)
        int_end = input_model.meta.exposure.integration_end  # one indexed
        if int_end is None:
            # Number of tables for the first (possibly only) spectral order.
            int_end = ntables_order[0]
            log.warning("INTEND not found; assuming a value of %d", int_end)

        # Columns of integration numbers & times of integration from the
        # INT_TIMES table.
        int_num = input_model.int_times["integration_number"]  # one indexed
        mid_utc = input_model.int_times["int_mid_MJD_UTC"]
        offset = int_start - int_num[0]
        if offset < 0 or int_end > int_num[-1]:
            log.warning(
                "Range of integration numbers in science data extends "
                "outside the range in INT_TIMES table."
            )
            log.warning("Can't use INT_TIMES table.")
            del int_num, mid_utc
            nrows = 0  # flag as bad
        else:
            log.debug("Times are from the INT_TIMES table.")
            time_arr = np.zeros(ntables, dtype=np.float64)
            j0 = 0
            for k in range(norders):
                ntables_current = ntables_order[k]
                j1 = j0 + ntables_current
                time_arr[j0:j1] = mid_utc[offset : offset + ntables_current]
                j0 += ntables_current

            int_times = Time(time_arr, format="mjd", scale="utc")

    if nrows == 0:
        log.debug("Times were computed from EXPSTART and TGROUP.")
        # Compute the delta time of each integration
        dt_arr = np.zeros(ntables, dtype=np.float64)
        dt = input_model.meta.exposure.group_time * (input_model.meta.exposure.ngroups + 1)
        j0 = 0
        for k in range(norders):
            ntables_current = ntables_order[k]
            j1 = j0 + ntables_current
            dt_arr[j0:j1] = np.arange(1, 1 + ntables_current) * dt - (dt / 2.0)
            j0 += ntables_current
        int_dt = TimeDelta(dt_arr, format="sec")

        # Compute the absolute time at the mid-point of each integration
        int_times = Time(input_model.meta.exposure.start_time, format="mjd") + int_dt

    # Store the times and flux sums in the table
    tbl["MJD"] = int_times.mjd
    tbl["whitelight_flux"] = fluxsums

    return tbl
