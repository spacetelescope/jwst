"""Sum the flux over all wavelengths in each integration as a function of time for the target."""

import logging

import numpy as np
from collections import OrderedDict
from astropy.table import QTable

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def white_light(input_model, min_wave=None, max_wave=None):
    """
    Compute the integrated flux over all wavelengths for a multi-integration extracted spectrum.

    Parameters
    ----------
    input_model : TSOMultiSpecModel
        Datamodel containing the multi-integration data.

    min_wave : float, optional
        Default wavelength minimum for integration.

    max_wave : float, optional
        Default wavelength maximum for integration.

    Returns
    -------
    tbl : astropy.table.table.QTable
        Table containing the integrated flux as a function of time.
    """
    if min_wave is None:
        min_wave = -1.0
    if max_wave is None:
        max_wave = 1.0e10

    # The input should contain one spectrum for each spectral
    # order.  NIRISS SOSS data can contain up to three orders.
    # Each row in the table is an integration.
    sporders = []  # list of spectral orders available
    order_list = []
    mid_times = None
    mid_tdbs = []
    flux_sums = None

    # Loop over the spectra in the input model and find mid times and fluxes
    for spec in input_model.spec:
        n_spec = len(spec.spec_table)

        # Figure out the spectral order for this spectrum
        spectral_order = getattr(spec, "spectral_order", None)
        sporders.append(spectral_order)
        order_list.extend([spectral_order] * n_spec)

        # Get mid times for all integrations in this order
        mid_time = spec.spec_table["MJD-AVG"]
        mid_tdb = spec.spec_table["TDB-MID"]
        good = ~np.isnan(mid_time)
        if mid_times is None:
            mid_times = mid_time
            mid_tdbs = mid_tdb
        else:
            mid_times = np.hstack([mid_times, mid_time])
            mid_tdbs = np.hstack([mid_tdbs, mid_tdb])

        # Create a wavelength mask, using cutoffs if specified, then
        # compute the flux sum for each integration in the input.
        wave_array = spec.spec_table["WAVELENGTH"]
        wave_mask = (wave_array >= min_wave) & (wave_array <= max_wave) & good[:, None]
        masked_flux = spec.spec_table["FLUX"].copy()
        masked_flux[~wave_mask] = np.nan
        flux_sum = np.nansum(masked_flux, axis=1)
        if flux_sums is None:
            flux_sums = flux_sum
        else:
            flux_sums = np.hstack([flux_sums, flux_sum])

        problems = np.sum(~good)
        if problems > 0:
            log.warning(
                f"There were {problems} spectra in order {spectral_order} "
                f"with no mid time ({100.0 * problems / n_spec} percent of spectra). "
                "These spectra will be ignored in the output table."
            )

    # Set up output table, removing problems
    tbl = _make_empty_output_table(input_model)
    good = ~np.isnan(mid_times)
    mid_times = mid_times[good]
    mid_tdbs = mid_tdbs[good]
    unique_mid_times, unq_indices = np.unique(mid_times, return_index=True)
    bjd_times = mid_tdbs[unq_indices]
    tbl["int_mid_MJD_UTC"] = unique_mid_times
    tbl["int_mid_BJD_TDB"] = bjd_times

    # Loop over the spectral orders and make separate table columns
    # for fluxes in each order, ensuring times are aligned even if one order
    # is not observed at a given time.
    for order in sporders:
        is_this_order = order_list == order
        time_is_in_this_order = np.isin(unique_mid_times, mid_times[is_this_order])

        # NaN-pad fluxes for times not represented in this order
        fluxes = np.full(len(unique_mid_times), np.nan)
        fluxes[time_is_in_this_order] = flux_sums[is_this_order]
        if len(sporders) > 1:
            colname = f"whitelight_flux_order_{order}"
        else:
            colname = "whitelight_flux"
        tbl[colname] = fluxes

    return tbl


def _make_empty_output_table(input_model):
    """
    Create an empty output table with the same metadata as the input model.

    Parameters
    ----------
    input_model : MultiSpecModel
        Datamodel containing the multi-integration data

    Returns
    -------
    astropy.table.table.QTable
        Empty table with the same metadata as the input model.
    """
    tbl_meta = OrderedDict()
    tbl_meta["instrument"] = input_model.meta.instrument.name
    tbl_meta["detector"] = input_model.meta.instrument.detector
    tbl_meta["exp_type"] = input_model.meta.exposure.type
    tbl_meta["subarray"] = input_model.meta.subarray.name
    tbl_meta["filter"] = input_model.meta.instrument.filter
    tbl_meta["pupil"] = input_model.meta.instrument.pupil
    tbl_meta["target_name"] = input_model.meta.target.catalog_name
    return QTable(meta=tbl_meta)
