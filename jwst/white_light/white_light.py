"""Sum the flux over all wavelengths in each integration as a function of time for the target."""

import logging
from collections import OrderedDict

import numpy as np
from astropy.table import QTable

log = logging.getLogger(__name__)

__all__ = ["white_light"]


def white_light(input_model, waverange_table=None, min_wave=None, max_wave=None):
    """
    Compute the integrated flux over all wavelengths for a multi-integration extracted spectrum.

    Parameters
    ----------
    input_model : TSOMultiSpecModel
        Datamodel containing the multi-integration data.
    waverange_table : astropy.table.Table, optional
        Wavelength range information from the reference file. These ranges will be
        superseded by user-specified values if provided.
    min_wave : float, optional
        User-specified wavelength minimum for integration.
    max_wave : float, optional
        User-specified wavelength maximum for integration.

    Returns
    -------
    tbl : astropy.table.table.QTable
        Table containing the integrated flux as a function of time.
    """
    # The input should contain separate spectra for each spectral
    # order or detector.  NIRISS SOSS data can contain up to three orders;
    # NIRSpec BOTS can contain up to two detectors.
    # Each row in the table is an integration.
    sporders = []  # list of spectral orders available
    detectors = []  # list of detectors available
    order_list = []
    detector_list = []
    mid_times = []
    mid_tdbs = []
    flux_sums = []

    # Loop over the spectra in the input model and find mid times and fluxes
    for spec in input_model.spec:
        n_spec = len(spec.spec_table)

        # Figure out the spectral order for this spectrum
        spectral_order = getattr(spec, "spectral_order", None)
        if spectral_order not in sporders:
            sporders.append(spectral_order)
        order_list.extend([spectral_order] * n_spec)

        # Do the same for the detector
        detector = getattr(spec, "detector", None)
        if detector not in detectors:
            detectors.append(detector)
        detector_list.extend([detector] * n_spec)

        # Determine wavelength range from either user-specified values or ref file
        min_wave, max_wave = _determine_wavelength_range(
            spectral_order,
            input_model.meta.instrument.filter,
            waverange_table=waverange_table,
            min_wave=min_wave,
            max_wave=max_wave,
        )

        # Get mid times for all integrations in this order
        mid_time = spec.spec_table["MJD-AVG"]
        mid_tdb = spec.spec_table["TDB-MID"]

        # Check for unique time stamps: keep only the first
        _, unq_idx = np.unique(mid_time, return_index=True)
        is_unique = np.full(mid_time.shape, False)
        is_unique[unq_idx] = True
        mid_time[~is_unique] = np.nan

        # Store time arrays
        good = ~np.isnan(mid_time)
        if len(mid_times) == 0:
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
        if len(flux_sums) == 0:
            flux_sums = flux_sum
        else:
            flux_sums = np.hstack([flux_sums, flux_sum])

        problems = np.sum(~good)
        if problems > 0:
            log.warning(
                f"There were {problems} spectra in order {spectral_order} "
                "with no mid time or duplicate mid time "
                f"({100.0 * problems / n_spec} percent of spectra). "
            )
            log.warning("These spectra will be ignored in the output table.")

    # Set up output table, removing problems
    tbl = _make_empty_output_table(input_model)
    good = ~np.isnan(mid_times)
    mid_times = mid_times[good]
    mid_tdbs = mid_tdbs[good]

    flux_sums = flux_sums[good]
    order_list = np.array(order_list)[good]
    detector_list = np.array(detector_list)[good]

    # Get time stamps for each detector - they generally have different values.
    max_rows = 0
    mjd_utc = {}
    bjd_tdb = {}
    for detector in detectors:
        is_this_detector = detector_list == detector
        mjd_det = mid_times[is_this_detector]
        unique_mid_times, unq_indices = np.unique(mjd_det, return_index=True)
        mjd_utc[detector] = unique_mid_times
        bjd_tdb[detector] = mid_tdbs[is_this_detector][unq_indices]
        if unique_mid_times.size > max_rows:
            max_rows = unique_mid_times.size

    # Loop over the detectors and spectral orders and make separate table columns
    # for times in each detector and fluxes in each order
    for detector in detectors:
        # Add the time column, with NaN-padding just in case the detector
        # timestamps do not quite align
        detector_rows = mjd_utc[detector].size
        mjd = np.full(max_rows, np.nan)
        mjd[:detector_rows] = mjd_utc[detector]
        bjd = np.full(max_rows, np.nan)
        bjd[:detector_rows] = bjd_tdb[detector]

        if len(detectors) > 1 or str(detector).upper() in ["NRS1", "NRS2"]:
            # add the detector to the column name if there are more than 1,
            # or the detectors are for NIRSpec
            detector_name = f"_{detector}"
        else:
            detector_name = ""
        tbl[f"MJD_UTC{detector_name}"] = mjd
        tbl[f"BJD_TDB{detector_name}"] = bjd

        for order in sporders:
            is_this_column = (order_list == order) & (detector_list == detector)
            time_is_in_this_column = np.isin(mjd, mid_times[is_this_column])

            # NaN-pad columns for times not represented in this column
            fluxes = np.full(max_rows, np.nan)
            fluxes[time_is_in_this_column] = flux_sums[is_this_column]

            colname = "whitelight_flux"
            if len(sporders) > 1:
                # add the spectral order to the column name if there are more than 1
                colname += f"_order_{order}"
            tbl[f"{colname}{detector_name}"] = fluxes

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
    tbl_meta["exp_type"] = input_model.meta.exposure.type
    tbl_meta["subarray"] = input_model.meta.subarray.name
    tbl_meta["filter"] = input_model.meta.instrument.filter
    tbl_meta["pupil"] = input_model.meta.instrument.pupil
    tbl_meta["target_name"] = input_model.meta.target.catalog_name
    return QTable(meta=tbl_meta)


def _determine_wavelength_range(order, filt, waverange_table=None, min_wave=None, max_wave=None):
    """
    Figure out wavelength range for a given filter and spectral order.

    If user-specified min/max wavelengths are provided, they will be used.
    Otherwise, the function will look for the wavelength range in the reference file info.
    If no reference file range is found, broad fallback values are returned,
    such that all wavelengths found in the input model spectral tables are used
    by the step.

    Parameters
    ----------
    order : int
        Spectral order for which to determine the wavelength range.
    filt : str
        Filter for which to determine the wavelength range.
    waverange_table : astropy.table.Table, optional
        Wavelength range information from the reference file.
    min_wave : float, optional
        User-specified minimum wavelength for integration.
    max_wave : float, optional
        User-specified maximum wavelength for integration.

    Returns
    -------
    tuple
        Minimum and maximum wavelengths for integration.
    """
    if waverange_table is not None:
        this_one = (waverange_table["order"] == int(order)) & (waverange_table["filter"] == filt)
        if not np.any(this_one):
            raise ValueError(
                f"No reference wavelength range found for order {order} and filter {filt}."
            )
        if np.sum(this_one) > 1:
            raise ValueError(
                f"Multiple reference wavelength ranges found for order {order} and filter {filt}."
            )
        min_wave_ref = waverange_table["min_wave"][this_one]
        max_wave_ref = waverange_table["max_wave"][this_one]
    else:
        # Set default values to well beyond any possible observed wavelength range
        min_wave_ref = -1.0
        max_wave_ref = 1.0e10

    # Only use ref file values if the user has not specified any values
    if min_wave is None:
        min_wave = min_wave_ref
    if max_wave is None:
        max_wave = max_wave_ref

    return min_wave, max_wave
