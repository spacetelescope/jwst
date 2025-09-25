import logging

import numpy as np
from astropy.nddata.bitmask import bitfield_to_boolean_mask
from scipy.interpolate import CubicSpline, UnivariateSpline
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import SossWaveGridModel, dqflags

from jwst.datamodels.utils.tso_multispec import make_tso_specmodel
from jwst.extract_1d.extract import populate_time_keywords
from jwst.extract_1d.soss_extract.atoca import ExtractionEngine, MaskOverlapError
from jwst.extract_1d.soss_extract.atoca_utils import (
    WebbKernel,
    get_wave_p_or_m,
    grid_from_map_with_extrapolation,
    make_combined_adaptive_grid,
    oversample_grid,
    throughput_soss,
)
from jwst.extract_1d.soss_extract.pastasoss import (
    CUTOFFS,
    _find_spectral_order_index,
    _verify_requested_orders,
    get_soss_wavemaps,
)
from jwst.extract_1d.soss_extract.soss_boxextract import (
    box_extract,
    estim_error_nearest_data,
    get_box_weights,
)
from jwst.extract_1d.soss_extract.soss_syscor import make_background_mask, soss_background
from jwst.lib import pipe_utils

log = logging.getLogger(__name__)

"""Shortest wavelength to consider for each spectral order"""
SHORT_CUTOFF = [None, 0.58, 0.63]

"""
Wavelengths short of which to consider Order 2 to be well-separated from order 1.
Two values are defined. The first is the cutoff longwave of which should be included
in the combined extraction; the second is the cutoff shortwave of which should be
included in the extraction of just that order. Some overlap is helpful to avoid
numerical issues.
"""
ORDER2_SEPARATION_CUTOFF = [0.77, 0.95]

ORDER_STR_TO_INT = {f"Order {order}": order for order in [1, 2, 3]}


__all__ = ["get_ref_file_args", "run_extract1d"]


def get_ref_file_args(ref_files, orders_requested=None):
    """
    Prepare the reference files for the extraction engine.

    Parameters
    ----------
    ref_files : dict
        A dictionary of the reference file DataModels, along with values
        for the subarray and pwcpos, i.e. the pupil wheel position.
    orders_requested : list or None, optional
        A list of the spectral orders requested for extraction.
        If None, all orders in the pastasoss reference file are used.

    Returns
    -------
    tuple
        The reference file args used with the extraction engine:
        (wavemaps, specprofiles, throughputs, kernels)
    """
    pastasoss_ref = ref_files["pastasoss"]
    specprofile_ref = ref_files["spec_profiles"]
    speckernel_ref = ref_files["speckernel"]
    n_pix = 2 * speckernel_ref.meta.halfwidth + 1
    speckernel_wv_range = [np.min(speckernel_ref.wavelengths), np.max(speckernel_ref.wavelengths)]

    refmodel_orders = [int(trace.spectral_order) for trace in pastasoss_ref.traces]
    if orders_requested is None:
        orders_requested = refmodel_orders
    else:
        orders_requested = _verify_requested_orders(orders_requested, refmodel_orders)

    wavemaps, spectraces = get_soss_wavemaps(
        ref_files["pwcpos"],
        refmodel=pastasoss_ref,
        subarray=ref_files["subarray"],
        padsize=None,
        spectraces=True,
        orders_requested=orders_requested,
    )

    # Collect spectral profiles, wavemaps, throughputs, kernels for all the orders
    spec_profiles = []
    throughputs = []
    kernels = []
    for order in orders_requested:
        order_idx = _find_spectral_order_index(pastasoss_ref, order)
        wavemap = wavemaps[order_idx]
        specprofile = specprofile_ref.profile[order_idx].data

        # apply padding to make specprofile and wavemap have same shape
        prof_shape0, prof_shape1 = specprofile.shape
        wavemap_shape0, wavemap_shape1 = wavemap.shape
        if prof_shape0 != wavemap_shape0:
            pad0 = (prof_shape0 - wavemap_shape0) // 2
            if pad0 > 0:
                specprofile = specprofile[pad0:-pad0, :]
            elif pad0 < 0:
                wavemap = wavemap[pad0:-pad0, :]
        if prof_shape1 != wavemap_shape1:
            pad1 = (prof_shape1 - wavemap_shape1) // 2
            if pad1 > 0:
                specprofile = specprofile[:, pad1:-pad1]
            elif pad1 < 0:
                wavemap = wavemap[:, pad1:-pad1]
        wavemaps[order_idx] = wavemap
        spec_profiles.append(specprofile)

        # make throughput interpolator
        thru = throughput_soss(
            pastasoss_ref.throughputs[order_idx].wavelength[:],
            pastasoss_ref.throughputs[order_idx].throughput[:],
        )
        throughputs.append(thru)

        # Build a kernel for this order
        wv_cent = np.zeros(wavemap.shape[1])

        # Get central wavelength as a function of columns
        col, _, wv = _get_trace_1d(spectraces, order)
        wv_cent[col] = wv

        # Set invalid values to zero
        idx_invalid = ~np.isfinite(wv_cent)
        wv_cent[idx_invalid] = 0.0

        kernel = WebbKernel(speckernel_ref.wavelengths, speckernel_ref.kernels, wv_cent, n_pix)
        valid_wavemap = (speckernel_wv_range[0] <= wavemap) & (wavemap <= speckernel_wv_range[1])
        wavemap = np.where(valid_wavemap, wavemap, 0.0)
        kernels.append(kernel)

    out = {
        "wavemaps": wavemaps,
        "spectraces": spectraces,
        "spec_profiles": spec_profiles,
        "throughputs": throughputs,
        "kernels": kernels,
        "subarray": ref_files["subarray"],
    }
    return out


def _get_trace_1d(spectraces, order):
    """
    Get the x, y, wavelength of the trace after applying the transform.

    Parameters
    ----------
    spectraces : list[np.ndarray]
        The list of spectral traces from the pastasoss reference file, one per spectral order.
    order : int
        The spectral order for which to return the trace parameters.

    Returns
    -------
    xtrace, ytrace, wavetrace : array[float]
        The x, y and wavelength of the trace.
    """
    order_idx = order - 1
    spectrace = spectraces[order_idx]
    xtrace = np.arange(CUTOFFS[order_idx])

    # CubicSpline requires monotonically increasing x arr
    if spectrace[0][0] - spectrace[0][1] > 0:
        spectrace = np.flip(spectrace, axis=1)

    trace_interp_y = CubicSpline(spectrace[0], spectrace[1])
    trace_interp_wave = CubicSpline(spectrace[0], spectrace[2])
    ytrace = trace_interp_y(xtrace)
    wavetrace = trace_interp_wave(xtrace)
    return xtrace, ytrace, wavetrace


def _estim_flux_first_order(
    scidata_bkg, scierr, scimask, wave_maps, spat_pros, thrpts, mask_trace_profile, threshold=1e-4
):
    """
    Roughly estimate the underlying flux of the target spectrum.

    This is done by simply masking out order 2 and retrieving the flux from order 1.

    Parameters
    ----------
    scidata_bkg : array
        A single background subtracted NIRISS SOSS detector image.
    scierr : array
        The uncertainties corresponding to the detector image.
    scimask : array
        Pixel mask to apply to the detector image.
    wave_maps : list[array]
        The wavelength maps for each spectral order.
    spat_pros : list[array]
        The spatial profiles for each spectral order.
    thrpts : list[func]
        The throughput functions for each spectral order.
    mask_trace_profile : array[bool]
        Mask determining the aperture used for extraction.
        Set to False where the pixel should be extracted.
    threshold : float, optional:
        The pixels with an aperture[order 2] > `threshold` are considered contaminated
        and will be masked. Default is 1e-4.

    Returns
    -------
    func
        A spline estimator that provides the underlying flux as a function of wavelength
    """
    # Define wavelength grid based on order 1 only (so first index)
    wave_grid = grid_from_map_with_extrapolation(wave_maps[0], spat_pros[0], n_os=1)

    # Mask parts contaminated by order 2 based on its spatial profile
    mask = (spat_pros[1] >= threshold) | mask_trace_profile[0]

    # Init extraction without convolution kernel (so extract the spectrum at order 1 resolution)
    args = {
        "wavemaps": [wave_maps[0]],
        "spec_profiles": [spat_pros[0]],
        "throughputs": [thrpts[0]],
        "kernels": [None],
    }
    engine = ExtractionEngine(args, wave_grid, [mask], global_mask=scimask, orders=[1])

    # Extract estimate
    spec_estimate = engine(scidata_bkg, scierr)

    # Interpolate
    idx = np.isfinite(spec_estimate)
    return UnivariateSpline(wave_grid[idx], spec_estimate[idx], k=3, s=0, ext=0)


def _get_native_grid_from_trace(spectraces, spectral_order):
    """
    Make a 1d-grid of the pixels boundary based on the wavelength solution.

    Parameters
    ----------
    spectraces : list[np.ndarray]
        The list of spectral traces from the pastasoss reference file, one per spectral order.
    spectral_order : int
        The spectral order for which to return the trace parameters.

    Returns
    -------
    wave : array[float]
        Grid of the pixels boundaries at the native sampling (1d array)
    col : array[int]
        The column number of the pixel
    """
    # From wavelength solution
    col, _, wave = _get_trace_1d(spectraces, spectral_order)

    # Keep only valid solution ...
    idx_valid = np.isfinite(wave)
    # ... and should correspond to subsequent columns
    is_subsequent = np.diff(col[idx_valid]) == 1
    if not is_subsequent.all():
        msg = f"Wavelength solution for order {spectral_order} contains gaps."
        log.warning(msg)
    wave = wave[idx_valid]
    col = col[idx_valid]
    log.debug(f"Wavelength range for order {spectral_order}: ({wave[[0, -1]]})")

    # Sort
    idx_sort = np.argsort(wave)
    wave = wave[idx_sort]
    col = col[idx_sort]

    return wave, col


def _get_grid_from_trace(spectraces, spectral_order, n_os):
    """
    Make a 1d-grid of the pixels boundary based on the wavelength solution.

    Parameters
    ----------
    spectraces : list[np.ndarray]
        The list of spectral traces from the pastasoss reference file, one per spectral order.
    spectral_order : int
        The spectral order for which to return the trace parameters.
    n_os : int or array
        The oversampling factor of the wavelength grid used when solving for
        the uncontaminated flux.

    Returns
    -------
    array[float]
        Grid of the pixels boundaries at the native sampling (1d array)
    """
    wave, _ = _get_native_grid_from_trace(spectraces, spectral_order)

    # Use pixel boundaries instead of the center values
    wv_upper_bnd, wv_lower_bnd = get_wave_p_or_m(wave[None, :])
    # `get_wave_p_or_m` returns 2d array, so keep only 1d
    wv_upper_bnd, wv_lower_bnd = wv_upper_bnd[0], wv_lower_bnd[0]
    # Each upper boundary should correspond the the lower boundary
    # of the following pixel, so only need to add the last upper boundary to complete the grid
    wv_upper_bnd, wv_lower_bnd = np.sort(wv_upper_bnd), np.sort(wv_lower_bnd)
    wave_grid = np.append(wv_lower_bnd, wv_upper_bnd[-1])

    # Oversample as needed
    return oversample_grid(wave_grid, n_os=n_os)


def _make_decontamination_grid(spectraces, rtol, max_grid_size, estimate, n_os):
    """
    Create the grid to use for the simultaneous extraction of order 1 and 2.

    The grid is made by:
    1) requiring that it satisfies the oversampling n_os
    2) trying to reach the specified tolerance for the spectral range shared between order 1 and 2
    3) trying to reach the specified tolerance in the rest of spectral range
    The max_grid_size overrules steps 2) and 3), so the precision may not be reached if
    the grid size needed is too large.

    Parameters
    ----------
    spectraces : list[np.ndarray]
        The list of spectral traces from the pastasoss reference file, one per spectral order.
    rtol : float
        The relative tolerance needed on a pixel model.
    max_grid_size : int
        Maximum grid size allowed.
    estimate : UnivariateSpline
        Estimate of the target flux as a function of wavelength in microns.
    n_os : int
        The oversampling factor of the wavelength grid used when solving for
        the uncontaminated flux.

    Returns
    -------
    wave_grid : 1d array
        The grid of the pixels boundaries at the native sampling.
    """
    # Build native grid for each  orders.
    spectral_orders = [2, 1]
    grids_ord = {}
    for sp_ord in spectral_orders:
        grids_ord[sp_ord] = _get_grid_from_trace(spectraces, sp_ord, n_os=n_os)

    # Build the list of grids given to make_combined_grid.
    # It must be ordered in increasing priority.
    # 1rst priority: shared wavelengths with order 1 and 2.
    # 2nd priority: remaining red part of order 1
    # 3rd priority: remaining blue part of order 2
    # So, split order 2 in 2 parts, the shared wavelength and the bluemost part
    is_shared = grids_ord[2] >= np.min(grids_ord[1])
    # Make sure order 1 is not more in the blue than order 2
    cond = grids_ord[1] > np.min(grids_ord[2][is_shared])
    grids_ord[1] = grids_ord[1][cond]
    # And make grid list
    all_grids = [grids_ord[2][is_shared], grids_ord[1], grids_ord[2][~is_shared]]

    # Cut order 2 at 0.77 (not smaller than that)
    # because there is no contamination there. Can be extracted afterward.
    # In the red, no cut.
    wv_range = [ORDER2_SEPARATION_CUTOFF[0], np.max(grids_ord[1])]

    # Finally, build the list of corresponding estimates.
    # The estimate for the overlapping part is the order 1 estimate.
    # There is no estimate yet for the blue part of order 2, so give a flat spectrum.
    def flat_fct(wv):
        return np.ones_like(wv)

    all_estimates = [estimate, estimate, flat_fct]

    # Generate the combined grid
    kwargs = {"rtol": rtol, "max_total_size": max_grid_size, "max_iter": 30}
    return make_combined_adaptive_grid(all_grids, all_estimates, wv_range, **kwargs)


def _append_tiktests(test_a, test_b):
    out = {}
    for key in test_a:
        out[key] = np.append(test_a[key], test_b[key], axis=0)

    return out


def _populate_tikho_attr(spec, tiktests, idx, sp_ord):
    spec.spectral_order = sp_ord
    spec.meta.soss_extract1d.type = "TEST"
    spec.meta.soss_extract1d.chi2 = tiktests["chi2"][idx]
    spec.meta.soss_extract1d.chi2_soft_l1 = tiktests["chi2_soft_l1"][idx]
    spec.meta.soss_extract1d.chi2_cauchy = tiktests["chi2_cauchy"][idx]
    spec.meta.soss_extract1d.reg = np.nansum(tiktests["reg"][idx] ** 2)
    spec.meta.soss_extract1d.factor = tiktests["factors"][idx]
    spec.int_num = 0


def _f_to_spec(f_order, grid_order, ref_file_args, pixel_grid, mask, sp_ord):
    """
    Bin the flux to the pixel grid and build a SpecModel.

    Parameters
    ----------
    f_order : np.array
        The solution f_k of the linear system.
    grid_order : np.array
        The wavelength grid of the solution, usually oversampled compared to the pixel grid.
    ref_file_args : dict
        The reference file arguments used by the ExtractionEngine.
    pixel_grid : np.array
        The pixel grid to which the flux should be binned.
    mask : np.array
        The mask of the pixels to be extracted.
    sp_ord : int
        The spectral order of the flux.

    Returns
    -------
    spec : SpecModel
        The SpecModel containing the extracted spectrum.
    """
    # Build 1d spectrum integrated over pixels
    pixel_grid = pixel_grid[np.newaxis, :]
    args = {
        "wavemaps": [pixel_grid],
        "spec_profiles": [np.ones_like(pixel_grid)],
        "throughputs": ref_file_args["throughputs"].copy(),
        "kernels": ref_file_args["kernels"].copy(),
    }
    model = ExtractionEngine(
        args,
        wave_grid=grid_order,
        mask_trace_profile=[mask[np.newaxis, :]],
        orders=[sp_ord],
    )
    f_binned = model.rebuild(f_order, fill_value=np.nan)

    pixel_grid = np.squeeze(pixel_grid)
    f_binned = np.squeeze(f_binned)

    # Remove Nans to save space
    is_valid = np.isfinite(f_binned)
    table_size = np.sum(is_valid)
    out_table = np.zeros(table_size, dtype=datamodels.SpecModel().spec_table.dtype)
    out_table["WAVELENGTH"] = pixel_grid[is_valid]
    out_table["FLUX"] = f_binned[is_valid]
    spec = datamodels.SpecModel(spec_table=out_table)
    spec.spectral_order = sp_ord

    return spec


def _build_tracemodel_order(i_order, engine, spectraces, f_k, mask):
    """
    Build the trace model for a specific spectral order.

    Parameters
    ----------
    i_order : int
        The spectral order index.
    engine : ExtractionEngine
        The engine used to model combined orders 1 and 2.
    spectraces : _type_
        _description_
    f_k : np.ndarray
        The extracted flux for the spectral order.
    mask : np.ndarray[bool]
        The global mask of pixels to be modeled. Bad pixels in the science data
        should remain unmasked.

    Returns
    -------
    tracemodel_ord : np.ndarray
        The modeled detector image for the spectral order.
    spec_ord : SpecModel
        The SpecModel containing the extracted spectrum for the spectral order.
    """
    # Pre-convolve the extracted flux (f_k) at the order's resolution
    # so that the convolution matrix must not be re-computed.
    flux_order = engine.kernels[i_order].dot(f_k)

    # Then must take the grid after convolution (smaller)
    grid_order = engine.wave_grid_c(i_order)

    # Keep only valid values to make sure there will be no Nans in the order model
    idx_valid = np.isfinite(flux_order)
    grid_order, flux_order = grid_order[idx_valid], flux_order[idx_valid]

    # Spectral order
    sp_ord = i_order + 1

    # Build model of the order
    # Give the identity kernel to the Engine (so no convolution)
    kernel = [np.array([1.0])]
    args_order = {
        "wavemaps": [engine.wave_map[i_order]],
        "spec_profiles": [engine.trace_profile[i_order]],
        "throughputs": [engine.throughput_orig[i_order]],
        "kernels": kernel,
    }
    model = ExtractionEngine(
        args_order, wave_grid=grid_order, mask_trace_profile=[mask], orders=[sp_ord]
    )

    # Project on detector and save in dictionary
    tracemodel_ord = model.rebuild(flux_order, fill_value=np.nan)

    # Build 1d spectrum integrated over pixels
    pixel_wave_grid, valid_cols = _get_native_grid_from_trace(spectraces, sp_ord)
    spec_ord = _f_to_spec(
        flux_order,
        grid_order,
        args_order,
        pixel_wave_grid,
        np.all(mask, axis=0)[valid_cols],
        sp_ord,
    )

    return tracemodel_ord, spec_ord


def _build_null_spec_table(wave_grid, order):
    """
    Build a SpecModel of entirely bad values.

    Parameters
    ----------
    wave_grid : np.array
        Input wavelengths
    order : int
        Spectral order

    Returns
    -------
    spec : SpecModel
        Null SpecModel. Flux values are NaN, DQ flags are 1,
        but note that DQ gets overwritten at end of run_extract1d
    """
    cut = SHORT_CUTOFF[order - 1]
    if cut is None:
        wave_grid_cut = wave_grid
    else:
        wave_grid_cut = wave_grid[wave_grid > cut]
    spec = datamodels.SpecModel()
    spec.spectral_order = order
    spec.meta.soss_extract1d.type = "OBSERVATION"
    spec.meta.soss_extract1d.factor = np.nan
    spec.spec_table = np.zeros((wave_grid_cut.size,), dtype=datamodels.SpecModel().spec_table.dtype)
    spec.spec_table["WAVELENGTH"] = wave_grid_cut
    spec.spec_table["FLUX"] = np.empty(wave_grid_cut.size) * np.nan
    spec.spec_table["DQ"] = np.ones(wave_grid_cut.size)
    spec.validate()

    return spec


def _do_tiktests(engine, scidata_bkg, scierr, estimate, spectraces, global_mask):
    # Find the tikhonov factor.
    # Initial pass 8 orders of magnitude with 10 grid points.
    guess_factor = engine.estimate_tikho_factors(estimate)
    log_guess = np.log10(guess_factor)
    factors = np.logspace(log_guess - 4, log_guess + 4, 10)
    all_tests = engine.get_tikho_tests(factors, scidata_bkg, scierr)
    tikfac = engine.best_tikho_factor(all_tests, fit_mode="all")

    # Refine across 4 orders of magnitude.
    tikfac = np.log10(tikfac)
    factors = np.logspace(tikfac - 2, tikfac + 2, 20)
    tiktests = engine.get_tikho_tests(factors, scidata_bkg, scierr)
    tikfac = engine.best_tikho_factor(tiktests, fit_mode="d_chi2")
    all_tests = _append_tiktests(all_tests, tiktests)

    # Save spectra in a list of SingleSpecModels for optional output
    spec_list = []
    for order in [1, 2]:
        for idx in range(len(all_tests["factors"])):
            f_k = all_tests["solution"][idx, :]
            _, spec_ord = _build_tracemodel_order(order - 1, engine, spectraces, f_k, global_mask)
            _populate_tikho_attr(spec_ord, all_tests, idx, order)
            spec_ord.meta.soss_extract1d.color_range = "RED"
            spec_list.append(spec_ord)
    return tikfac, spec_list


def _compute_box_weights(spectraces, shape, width, orders_requested):
    """
    Determine the weights for the box extraction.

    Parameters
    ----------
    spectraces : list[np.ndarray]
        The list of spectral traces from the pastasoss reference file, one per spectral order.
    shape : tuple
        The shape of the detector image.
    width : int
        The width of the box aperture.
    orders_requested : list
        List of orders to be extracted.

    Returns
    -------
    box_weights : dict
        A dictionary of the weights for each order.
    wavelengths : dict
        A dictionary of the wavelengths for each order.
    """
    # Extract each order from order list
    box_weights, wavelengths = {}, {}
    order_str = {order: f"Order {order}" for order in orders_requested}
    for order_integer in orders_requested:
        # Order string-name is used more often than integer-name
        order = order_str[order_integer]

        log.debug(f"Compute box weights for {order}.")

        # Define the box aperture
        xtrace, ytrace, wavelengths[order] = _get_trace_1d(spectraces, order_integer)
        box_weights[order] = get_box_weights(ytrace, width, shape, cols=xtrace)

    return box_weights, wavelengths


def _model_single_order(
    data_order,
    err_order,
    ref_file_args,
    mask_fit,
    mask_rebuild,
    order,
    wave_grid,
    valid_cols,
    tikfac_log_range,
    save_tiktests=False,
):
    """
    Extract an output spectrum for a single spectral order using the ATOCA algorithm.

    The Tikhonov factor is derived in two stages: first, ten factors are tested
    spanning tikfac_log_range, and then a further 20 factors are tested across
    2 orders of magnitude in each direction around the best factor from the first stage.
    The best-fitting model and spectrum are reconstructed using the best-fit Tikhonov factor
    and respecting mask_rebuild.

    Parameters
    ----------
    data_order : np.array
        The 2D data array for the spectral order to be extracted.
    err_order : np.array
        The 2D error array for the spectral order to be extracted.
    ref_file_args : dict
        The reference file arguments used by the ExtractionEngine.
    mask_fit : np.array
        Mask determining the aperture used for extraction. This typically includes
        detector bad pixels and any pixels that are not part of the trace
    mask_rebuild : np.array
        Mask determining the aperture used for rebuilding the trace. This typically includes
        only pixels that do not belong to either spectral trace, i.e., regions of the detector
        where no real data could exist.
    order : int
        The spectral order to be extracted.
    wave_grid : np.array
        The wavelength grid used to model the data.
    valid_cols : np.array
        The columns of the detector that are valid for extraction.
    tikfac_log_range : list
        The range of Tikhonov factors to test, in log space.
    save_tiktests : bool, optional
        If True, save the intermediate models and spectra for each Tikhonov factor tested.

    Returns
    -------
    model : np.array
        Model derived from the best Tikhonov factor, same shape as data_order.
    spec_list : list of SpecModel
        If save_tiktests is True, returns a list of the model spectra
        for each Tikhonov factor tested,
        with the best-fitting spectrum last in the list.
        If save_tiktests is False, returns a one-element list with the best-fitting spectrum.

    Notes
    -----
    The last spectrum in the list of SpecModels lacks
    the "chi2", "chi2_soft_l1", "chi2_cauchy", and "reg" attributes,
    as these are only calculated for the intermediate models. The last spectrum is not
    necessarily identical to any of the spectra in the list, as it is reconstructed according to
    mask_rebuild instead of fit respecting mask_fit; that is, bad pixels are included.
    """

    # The throughput and kernel is not needed here
    # set them so they have no effect on the extraction.
    def throughput(wavelength):
        return np.ones_like(wavelength)

    kernel = np.array([1.0])
    ref_file_args["throughputs"] = [throughput]
    ref_file_args["kernels"] = [kernel]

    # Define wavelength grid with oversampling of 3 (should be enough)
    wave_grid_os = oversample_grid(wave_grid, n_os=3)

    # Initialize the Engine.
    engine = ExtractionEngine(
        ref_file_args,
        wave_grid=wave_grid_os,
        mask_trace_profile=[mask_fit],
        orders=[order],
    )

    # Find the tikhonov factor.
    # Initial pass with tikfac_range.
    factors = np.logspace(tikfac_log_range[0], tikfac_log_range[-1], 10)
    all_tests = engine.get_tikho_tests(factors, data_order, err_order)
    tikfac = engine.best_tikho_factor(tests=all_tests, fit_mode="all")

    # Refine across 4 orders of magnitude.
    tikfac = np.log10(tikfac)
    factors = np.logspace(tikfac - 2, tikfac + 2, 20)
    tiktests = engine.get_tikho_tests(factors, data_order, err_order)
    tikfac = engine.best_tikho_factor(tiktests, fit_mode="d_chi2")
    all_tests = _append_tiktests(all_tests, tiktests)

    # Run the extract method of the Engine.
    f_k_final = engine(data_order, err_order, tikhonov=True, factor=tikfac)

    # Save binned spectra in a list of SingleSpecModels for optional output
    spec_list = []
    if save_tiktests:
        for idx in range(len(all_tests["factors"])):
            f_k = all_tests["solution"][idx, :]

            # Build 1d spectrum integrated over pixels
            spec_ord = _f_to_spec(
                f_k,
                wave_grid_os,
                ref_file_args,
                wave_grid,
                np.all(mask_rebuild, axis=0)[valid_cols],
                order,
            )
            _populate_tikho_attr(spec_ord, all_tests, idx, order)

            # Add the result to spec_list
            spec_list.append(spec_ord)

    # Rebuild trace, including bad pixels
    engine = ExtractionEngine(
        ref_file_args,
        wave_grid=wave_grid_os,
        mask_trace_profile=[mask_rebuild],
        orders=[order],
    )
    model = engine.rebuild(f_k_final, fill_value=np.nan)

    # Build 1d spectrum integrated over pixels
    spec_ord = _f_to_spec(
        f_k_final,
        wave_grid_os,
        ref_file_args,
        wave_grid,
        np.all(mask_rebuild, axis=0)[valid_cols],
        order,
    )
    spec_ord.meta.soss_extract1d.factor = tikfac
    spec_ord.meta.soss_extract1d.type = "OBSERVATION"

    # Add the result to spec_list
    spec_list.append(spec_ord)
    return model, spec_list


class Integration:
    def __init__(
        self,
        scidata,
        scierr,
        scimask,
        refmask,
        ref_file_args,
        box_weights,
        do_bkgsub=True,
        extract_order3=True,
    ):
        self.scidata = scidata
        self.scierr = scierr
        self.scimask = scimask
        self.refmask = refmask
        self.box_weights = box_weights

        if extract_order3:
            self.order_list = [1, 2, 3]
        else:
            self.order_list = [1, 2]
        self.order_strs = [f"Order {order}" for order in self.order_list]
        self.order_indices = [o - 1 for o in self.order_list]

        self._validate_masks()
        self._subtract_bkg(do_bkgsub)

        # unpack ref file args
        self.wavemaps = ref_file_args["wavemaps"]
        self.spectraces = ref_file_args["spectraces"]
        self.spec_profiles = ref_file_args["spec_profiles"]
        self.throughputs = ref_file_args["throughputs"]
        self.kernels = ref_file_args["kernels"]
        self.subarray = ref_file_args["subarray"]

    def _subtract_bkg(self, do_bkgsub):
        # Perform background correction.
        if do_bkgsub:
            log.info("Applying background subtraction.")
            bkg_mask = make_background_mask(self.scidata, width=40)
            self.scidata_bkg, self.col_bkg = soss_background(self.scidata, self.scimask, bkg_mask)
        else:
            log.info("Skip background subtraction.")
            self.scidata_bkg = self.scidata
            self.col_bkg = np.zeros(self.scidata.shape[1])

    def _validate_masks(self):
        # Make sure there aren't any nans not flagged in scimask
        not_finite = ~(np.isfinite(self.scidata) & np.isfinite(self.scierr))
        if (not_finite & ~self.scimask).any():
            log.warning(
                "Input contains invalid values that "
                "are not flagged correctly in the dq map. "
                "Masking those values for extraction."
            )
            self.scimask |= not_finite
            self.refmask &= ~not_finite

    def model_image(
        self,
        wave_grid=None,
        estimate=None,
        tikfac=None,
        rtol=1e-4,
        n_os=2,
        max_grid_size=20000,
        threshold=1e-2,
    ):
        # Some error values are 0, we need to mask those pixels for the extraction engine.
        scimask = self.scimask | ~(self.scierr > 0)

        # Define mask based on box aperture
        # (we want to model each contaminated pixels that will be extracted)
        mask_trace_profile = [
            (~(self.box_weights[self.order_strs[i]] > 0)) | (self.refmask)
            for i in self.order_indices
        ]

        # Define mask of pixel to model (all pixels inside box aperture)
        global_mask = np.all(mask_trace_profile, axis=0).astype(bool)

        # Rough estimate of the underlying flux
        if (tikfac is None or wave_grid is None) and estimate is None:
            estimate = _estim_flux_first_order(
                self.scidata_bkg,
                self.scierr,
                scimask,
                self.wavemaps,
                self.spec_profiles,
                self.throughputs,
                mask_trace_profile,
            )

        # Generate grid based on estimate if not given
        if wave_grid is None:
            log.info(f"wave_grid not given: generating grid based on rtol={rtol}")
            wave_grid = _make_decontamination_grid(
                self.spectraces, rtol, max_grid_size, estimate, n_os
            )
            log.debug(
                f"wave_grid covering from {wave_grid.min()} to {wave_grid.max()}"
                f" with {wave_grid.size} points"
            )
        else:
            log.info("Using previously computed or user specified wavelength grid.")

        # Initialize the Engine for combined extraction of orders 1 and 2
        args_orders_12 = {
            "wavemaps": self.wavemaps[:2],
            "spec_profiles": self.spec_profiles[:2],
            "throughputs": self.throughputs[:2],
            "kernels": self.kernels[:2],
        }
        engine = ExtractionEngine(
            args_orders_12,
            wave_grid=wave_grid,
            mask_trace_profile=mask_trace_profile[:2],
            global_mask=scimask,
            threshold=threshold,
            orders=[1, 2],
        )

        if tikfac is None:
            log.info("Solving for the optimal Tikhonov factor.")
            save_tiktests = True
            tikfac, spec_list = _do_tiktests(
                engine,
                self.scidata_bkg,
                self.scierr,
                estimate,
                self.spectraces,
                global_mask,
            )

        else:
            save_tiktests = False
            spec_list = []

        log.info(f"Using a Tikhonov factor of {tikfac}")

        # Run the extract method of the Engine.
        f_k = engine(self.scidata_bkg, self.scierr, tikhonov=True, factor=tikfac)

        # Compute the log-likelihood of the best fit.
        logl = engine.compute_likelihood(f_k, self.scidata_bkg, self.scierr)

        log.info(f"Optimal solution has a log-likelihood of {logl}")

        # Create a new instance of the engine for evaluating the trace model.
        # This allows bad pixels and pixels below the threshold to be reconstructed as well.
        # Model the traces for each order separately.
        tracemodels = {}
        for i_order, order in enumerate([1, 2]):
            log.debug(f"Building the model image of {order}.")

            tracemodel_ord, spec_ord = _build_tracemodel_order(
                i_order, engine, self.spectraces, f_k, global_mask
            )
            spec_ord.meta.soss_extract1d.factor = tikfac
            spec_ord.meta.soss_extract1d.color_range = "RED"
            spec_ord.meta.soss_extract1d.type = "OBSERVATION"

            # Project on detector and save in dictionary
            tracemodels[self.order_strs[i_order]] = tracemodel_ord

            # Add the result to spec_list
            spec_list.append(spec_ord)

        # Make a null tracemodel to be overwritten later for order 3
        if 3 in self.order_list:
            tracemodels["Order 3"] = np.zeros_like(tracemodels["Order 2"]) * np.nan

        # Model the blue part of order 2 and all of order 3 assuming they are well-separated
        # from order 1
        for order in self.order_list[1:]:
            if self.subarray == "SUBSTRIP96":
                continue
            if order not in self.order_list:
                continue
            if order == 2:
                cutoff = ORDER2_SEPARATION_CUTOFF[1]
            else:
                cutoff = None

            idx_order = np.array(self.order_indices)[np.array(self.order_list) == order][0]
            order_str = self.order_strs[idx_order]
            log.info(f"Generate model for well-separated part of {order_str}")

            # Take only the second order's specific ref_files
            ref_file_order = {
                "wavemaps": [self.wavemaps[idx_order]],
                "spec_profiles": [self.spec_profiles[idx_order]],
                "throughputs": [self.throughputs[idx_order]],
                "kernels": [self.kernels[idx_order]],
            }

            # Mask for the fit. All valid pixels inside box aperture
            mask_fit = mask_trace_profile[idx_order] | scimask

            # Build 1d spectrum integrated over pixels
            pixel_wave_grid, valid_cols = _get_native_grid_from_trace(self.spectraces, order)

            # Hardcode wavelength highest boundary as well.
            # Must overlap with lower limit in make_decontamination_grid
            if cutoff is not None:
                is_in_wv_range = pixel_wave_grid < cutoff
                pixel_wave_grid, valid_cols = (
                    pixel_wave_grid[is_in_wv_range],
                    valid_cols[is_in_wv_range],
                )

            # Range of initial tikhonov factors
            tikfac_log_range = np.log10(tikfac) + np.array([-2, 8])

            # Model with atoca
            try:
                model, spec_ord = _model_single_order(
                    self.scidata_bkg,
                    self.scierr,
                    ref_file_order,
                    mask_fit,
                    global_mask,
                    order,
                    pixel_wave_grid,
                    valid_cols,
                    tikfac_log_range,
                    save_tiktests=save_tiktests,
                )

            except MaskOverlapError:
                log.error(
                    "Not enough unmasked pixels to model the remaining part of order 2."
                    " Model and spectrum will be NaN in that spectral region."
                )
                spec_ord = [_build_null_spec_table(pixel_wave_grid, order)]
                model = np.nan * np.ones_like(self.scidata_bkg)

            # Keep only pixels from which order 2 contribution
            # is not already modeled.
            already_modeled = np.isfinite(tracemodels[order_str])
            model = np.where(already_modeled, 0.0, model)

            # Add to tracemodels
            both_nan = np.isnan(tracemodels[order_str]) & np.isnan(model)
            tracemodels[order_str] = np.nansum([tracemodels[order_str], model], axis=0)
            tracemodels[order_str][both_nan] = np.nan

            # Add the result to spec_list
            for sp in spec_ord:
                if order == 2:
                    sp.meta.soss_extract1d.color_range = "BLUE"
                else:
                    sp.meta.soss_extract1d.color_range = "ALL"
            spec_list += spec_ord

        return tracemodels, tikfac, logl, wave_grid, spec_list

    def decontaminate_image(self, tracemodels):
        """
        Perform decontamination of the image based on the trace models.

        Parameters
        ----------
        tracemodels : dict
            Dictionary of the modeled detector images for each order.

        Returns
        -------
        decontaminated_data : dict
            Dictionary of the decontaminated data for each order.
        """
        log.debug("Performing the decontamination.")
        # List of modeled orders
        mod_order_list = tracemodels.keys()

        # Extract each order from order list
        decontaminated_data = {}
        for order in self.order_strs:
            # Decontaminate using all other modeled orders
            decont = self.scidata_bkg.copy()
            for mod_order in mod_order_list:
                if mod_order != order:
                    log.debug(f"Decontaminating {order} from {mod_order} using model.")
                    is_valid = np.isfinite(tracemodels[mod_order])
                    decont = decont - np.where(is_valid, tracemodels[mod_order], 0.0)

            # Save results
            decontaminated_data[order] = decont
        return decontaminated_data

    def extract_image(self, decontaminated_data, bad_pix="model", tracemodels=None):
        """
        Perform the box-extraction on the image using the trace model to correct for contamination.

        Parameters
        ----------
        decontaminated_data : array[float]
            A single background subtracted NIRISS SOSS detector image.
        bad_pix : str
            How to handle the bad pixels. Options are 'masking' and 'model'.
            'masking' will simply mask the bad pixels, such that the number of pixels
            in each column in the box extraction will not be constant, while the
            'model' option uses `tracemodels` to replace the bad pixels.
        tracemodels : dict
            Dictionary of the modeled detector images for each order.

        Returns
        -------
        fluxes, fluxerrs, npixels : dict
            Each output is a dictionary, with each extracted order as a key.
        """
        # Init models with an empty dictionary if not given
        if tracemodels is None:
            tracemodels = {}

        # Which orders to extract (extract the ones with given box aperture).
        order_list = self.box_weights.keys()

        # Create dictionaries for the output spectra.
        fluxes, fluxerrs, npixels = {}, {}, {}

        log.info("Performing the box extraction.")

        # Extract each order from order list
        for order in order_list:
            box_w_ord = self.box_weights[order]
            decont = decontaminated_data[order]
            # Replace bad pixels with trace model
            if (bad_pix == "model") and (order in list(tracemodels.keys())):
                # Some pixels might not be modeled by the bad pixel models
                is_modeled = np.isfinite(tracemodels[order])
                # Replace bad pixels
                decont = np.where(self.scimask & is_modeled, tracemodels[order], decont)

                # Replace error estimate of the bad pixels
                # using other valid pixels of similar value.
                # The pixel to be estimated are the masked pixels in the region of extraction
                # with available model.
                extraction_region = box_w_ord > 0
                pix_to_estim = extraction_region & self.scimask & is_modeled
                # Use only valid pixels (not masked) in the extraction region
                # for the empirical estimation
                valid_pix = extraction_region & ~self.scimask
                scierr_ord = estim_error_nearest_data(self.scierr, decont, pix_to_estim, valid_pix)

                # Update the scimask for box extraction:
                # the pixels that are modeled are not masked anymore, so set to False.
                # Note that they have to be in the extraction region
                # to ensure that scierr is also valid
                scimask_ord = np.where(is_modeled, False, self.scimask)
                log.info(f"Bad pixels in {order} are replaced with trace model.")

            else:
                scimask_ord = self.scimask
                scierr_ord = self.scierr
                log.info(
                    f"Bad pixels in {order} will be masked instead of modeled: "
                    "Trace model unavailable or not requested."
                )

            # Perform the box extraction and save
            out = box_extract(decont, scierr_ord, scimask_ord, box_w_ord)
            _, fluxes[order], fluxerrs[order], npixels[order] = out

        return fluxes, fluxerrs, npixels


def _process_one_integration(
    scidata,
    scierr,
    scimask,
    refmask,
    ref_file_args,
    box_weights,
    wavelengths,
    soss_kwargs,
    wave_grid=None,
    tikfac=None,
    generate_model=True,
    int_num=None,
):
    integration = Integration(
        scidata,
        scierr,
        scimask,
        refmask,
        ref_file_args,
        box_weights,
        extract_order3=soss_kwargs["order_3"],
        do_bkgsub=soss_kwargs["subtract_background"],
    )

    # Model the traces based on optics filter configuration (CLEAR or F277W)
    if generate_model:
        # Model the image.
        tracemodels, tikfac, _, wave_grid, atoca_list = integration.model_image(
            estimate=soss_kwargs["estimate"],
            wave_grid=wave_grid,
            tikfac=tikfac,
            rtol=soss_kwargs["rtol"],
            n_os=soss_kwargs["n_os"],
            max_grid_size=soss_kwargs["max_grid_size"],
            threshold=soss_kwargs["threshold"],
        )

        # Add atoca spectra to multispec for output
        for i, spec in enumerate(atoca_list):
            # If it was a test, not the best spectrum,
            # int_num is already set to 0.
            if not spec.hasattr("int_num") and int_num is not None:
                spec.int_num = int_num
            atoca_list[i] = spec
    else:
        # Return empty tracemodels
        tracemodels = {}
        atoca_list = []

    # Decontaminate the data using trace models (if tracemodels not empty)
    data_to_extract = integration.decontaminate_image(tracemodels)

    # Use the bad pixel models to perform a de-contaminated extraction.
    result = integration.extract_image(
        data_to_extract,
        bad_pix=soss_kwargs["bad_pix"],
        tracemodels=tracemodels,
    )
    fluxes, fluxerrs, npixels = result

    # Save trace models for output reference
    for order in tracemodels:
        # Put NaNs to zero
        model_ord = tracemodels[order]
        model_ord = np.where(np.isfinite(model_ord), model_ord, 0.0)
        tracemodels[order] = model_ord

    # Copy spectral data for each order into the output model.
    spec_list = {}
    for order in fluxes.keys():
        table_size = len(wavelengths[order])

        out_table = np.zeros(table_size, dtype=datamodels.SpecModel().spec_table.dtype)
        out_table["WAVELENGTH"] = wavelengths[order][:table_size]
        out_table["FLUX"] = fluxes[order][:table_size]
        out_table["FLUX_ERROR"] = fluxerrs[order][:table_size]
        out_table["DQ"] = np.zeros(table_size)
        out_table["BACKGROUND"] = integration.col_bkg[:table_size]
        out_table["NPIXELS"] = npixels[order][:table_size]

        spec = datamodels.SpecModel(spec_table=out_table)

        # Add integration number and spectral order
        spec.spectral_order = ORDER_STR_TO_INT[order]
        if int_num is not None:
            spec.int_num = int_num

        spec_list[order] = spec

    return tracemodels, spec_list, atoca_list, tikfac, wave_grid


def run_extract1d(
    input_model,
    pastasoss_ref_name,
    specprofile_ref_name,
    speckernel_ref_name,
    subarray,
    soss_filter,
    soss_kwargs,
):
    """
    Run the spectral extraction on NIRISS SOSS data.

    Parameters
    ----------
    input_model : DataModel
        The input DataModel.
    pastasoss_ref_name : str
        Name of the pastasoss reference file.
    specprofile_ref_name : str
        Name of the specprofile reference file.
    speckernel_ref_name : str
        Name of the speckernel reference file.
    subarray : str
        Subarray on which the data were recorded; one of 'SUBSTRIP96',
        'SUBSTRIP256' or 'FULL'.
    soss_filter : str
        Filter in place during observations; one of 'CLEAR' or 'F277W'.
    soss_kwargs : dict
        Dictionary of keyword arguments passed from extract_1d_step.

    Returns
    -------
    output_model : DataModel
        DataModel containing the extracted spectra.
    """
    # Generate the atoca models or not (not necessarily for decontamination)
    generate_model = soss_kwargs["atoca"] or (soss_kwargs["bad_pix"] == "model")

    if soss_filter != "CLEAR" and generate_model:
        # No model can be fit for F277W yet, missing throughput reference files.
        msg = f"No extraction possible for filter {soss_filter}."
        log.critical(msg)
        raise ValueError(msg)

    # Read the reference files.
    pastasoss_ref = datamodels.PastasossModel(pastasoss_ref_name)
    specprofile_ref = datamodels.SpecProfileModel(specprofile_ref_name)
    speckernel_ref = datamodels.SpecKernelModel(speckernel_ref_name)

    # Map the order integer names to the string names
    if (soss_kwargs["order_3"]) and (subarray != "SUBSTRIP96"):
        order_list = [1, 2, 3]
    else:
        order_list = [1, 2]
    refmodel_orders = [int(trace.spectral_order) for trace in pastasoss_ref.traces]
    order_list = _verify_requested_orders(order_list, refmodel_orders)

    ref_files = {}
    ref_files["pastasoss"] = pastasoss_ref
    ref_files["spec_profiles"] = specprofile_ref
    ref_files["speckernel"] = speckernel_ref
    ref_files["subarray"] = subarray
    ref_files["pwcpos"] = input_model.meta.instrument.pupil_position

    # Unpack wave_grid if wave_grid_in was specified.
    wave_grid_in = soss_kwargs["wave_grid_in"]
    if wave_grid_in is not None:
        log.info(f"Loading wavelength grid from {wave_grid_in}.")
        wave_grid = datamodels.SossWaveGridModel(wave_grid_in).wavegrid
        # Make sure it has the correct precision
        wave_grid = wave_grid.astype("float64")
    else:
        # wave_grid will be estimated later in the first call of `_model_image`
        log.info("Wavelength grid was not specified. Setting `wave_grid` to None.")
        wave_grid = None

    # Convert estimate to cubic spline if given.
    # It should be a SpecModel or a file name (string)
    estimate = soss_kwargs.pop("estimate")
    if estimate is not None:
        log.info("Converting the estimate of the flux to spline function.")

        # Convert estimate to cubic spline
        estimate = datamodels.open(estimate)
        wv_estimate = estimate.spec_table["WAVELENGTH"]
        flux_estimate = estimate.spec_table["FLUX"]
        # Keep only finite values
        idx = np.isfinite(flux_estimate)
        estimate = UnivariateSpline(wv_estimate[idx], flux_estimate[idx], k=3, s=0, ext=0)
    soss_kwargs["estimate"] = estimate

    # Initialize the output model.
    output_model = datamodels.TSOMultiSpecModel()
    output_model.update(input_model)  # Copy meta data from input to output.

    # Initialize output spectra returned by ATOCA
    # NOTE: these diagnostic spectra are formatted as a simple MultiSpecModel,
    # with integrations in separate spectral extensions.
    output_atoca = datamodels.MultiSpecModel()
    output_atoca.update(input_model)

    # Initialize output references (model of the detector and box aperture weights).
    output_references = datamodels.SossExtractModel()
    output_references.update(input_model)

    # Convert to Cube if datamodels is an ImageModel
    if isinstance(input_model, datamodels.ImageModel):
        cube_model = datamodels.CubeModel(shape=(1, *input_model.shape))
        cube_model.data = input_model.data[None, :, :]
        cube_model.err = input_model.err[None, :, :]
        cube_model.dq = input_model.dq[None, :, :]
        nimages = 1
        log.info("Input is an ImageModel, processing a single integration.")
    elif isinstance(input_model, datamodels.CubeModel):
        cube_model = input_model
        nimages = len(cube_model.data)
        log.info(f"Input is a CubeModel containing {nimages} integrations.")
    else:
        msg = "Only ImageModel and CubeModel are implemented for the NIRISS SOSS extraction."
        log.critical(msg)
        raise TypeError(msg)

    # Prepare the reference file arguments.
    ref_file_args = get_ref_file_args(ref_files, orders_requested=order_list)

    # Run the first integration to get the Tikhonov factor for the rest
    scidata = cube_model.data[0].astype("float64")
    scierr = cube_model.err[0].astype("float64")
    scimask = np.bitwise_and(cube_model.dq[0], dqflags.pixel["DO_NOT_USE"]).astype(bool)
    refmask = bitfield_to_boolean_mask(
        cube_model.dq[0], ignore_flags=dqflags.pixel["REFERENCE_PIXEL"], flip_bits=True
    )

    # Pre-compute the weights for box extraction (used in modeling and extraction)
    box_weights, wavelengths = _compute_box_weights(
        ref_file_args["spectraces"],
        scidata.shape,
        width=soss_kwargs["width"],
        orders_requested=order_list,
    )
    # FIXME: hardcoding the substrip96 weights to unity is a band-aid solution
    if subarray == "SUBSTRIP96":
        box_weights["Order 2"] = np.ones((96, 2048))

    tracemodels, spec_list, atoca_list, tikfac_first, wave_grid_first = _process_one_integration(
        scidata,
        scierr,
        scimask,
        refmask,
        ref_file_args,
        box_weights,
        wavelengths,
        soss_kwargs,
        wave_grid=wave_grid,
        tikfac=soss_kwargs["tikfac"],
        generate_model=generate_model,
        int_num=0,
    )
    for atoca_spec in atoca_list:
        output_atoca.spec.append(atoca_spec)

    # Loop over images.
    all_tracemodels = {order: [tracemodels[order]] for order in tracemodels}
    output_spec_list = {order: [spec_list[order]] for order in spec_list}
    for i in range(1, nimages):
        log.info(f"Processing integration {i + 1} of {nimages}.")

        # Unpack the i-th image, set dtype to float64 and convert DQ to boolean mask.
        scidata = cube_model.data[i].astype("float64")
        scierr = cube_model.err[i].astype("float64")
        scimask = np.bitwise_and(cube_model.dq[i], dqflags.pixel["DO_NOT_USE"]).astype(bool)
        refmask = bitfield_to_boolean_mask(
            cube_model.dq[i], ignore_flags=dqflags.pixel["REFERENCE_PIXEL"], flip_bits=True
        )

        tracemodels, spec_list, atoca_list, _, _ = _process_one_integration(
            scidata,
            scierr,
            scimask,
            refmask,
            ref_file_args,
            box_weights,
            wavelengths,
            soss_kwargs,
            wave_grid=wave_grid_first,
            tikfac=tikfac_first,
            generate_model=generate_model,
            int_num=i,
        )
        for order in tracemodels:
            all_tracemodels[order].append(tracemodels[order])
        for order in spec_list:
            output_spec_list[order].append(spec_list[order])

    # Make a TSOSpecModel from the output spec list
    for order in output_spec_list:
        tso_spec = make_tso_specmodel(
            output_spec_list[order], segment=input_model.meta.exposure.segment_number
        )
        output_model.spec.append(tso_spec)

    # Update output model
    output_model.meta.soss_extract1d.width = soss_kwargs["width"]
    output_model.meta.soss_extract1d.apply_decontamination = soss_kwargs["atoca"]
    output_model.meta.soss_extract1d.tikhonov_factor = soss_kwargs["tikfac"]
    output_model.meta.soss_extract1d.oversampling = soss_kwargs["n_os"]
    output_model.meta.soss_extract1d.threshold = soss_kwargs["threshold"]
    output_model.meta.soss_extract1d.bad_pix = soss_kwargs["bad_pix"]

    # Save output references
    for order in all_tracemodels:
        # Convert from list to array
        tracemod_ord = np.array(all_tracemodels[order])
        # Save
        order_int = ORDER_STR_TO_INT[order]
        setattr(output_references, f"order{order_int}", tracemod_ord)

    for order in box_weights:
        # Convert from list to array
        box_w_ord = np.array(box_weights[order])
        # repeat along axis zero to have shape (nints, y, x)
        box_w_ord = np.repeat(box_w_ord[None, :, :], nimages, axis=0)
        # Save
        order_int = ORDER_STR_TO_INT[order]
        setattr(output_references, f"aperture{order_int}", box_w_ord)

    if pipe_utils.is_tso(input_model):
        log.info("Populating INT_TIMES keywords from input table.")
        populate_time_keywords(input_model, output_model)
        output_model.int_times = input_model.int_times.copy()

    if soss_kwargs["wave_grid_out"] is not None:
        wave_grid_model = SossWaveGridModel(wavegrid=wave_grid)
        log.info(f"Saving soss_wave_grid to {soss_kwargs['wave_grid_out']}")
        wave_grid_model.save(path=soss_kwargs["wave_grid_out"])
        wave_grid_model.close()

    return output_model, output_references, output_atoca
