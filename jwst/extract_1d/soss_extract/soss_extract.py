import logging
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from astropy.nddata.bitmask import bitfield_to_boolean_mask
from astropy.utils.decorators import lazyproperty
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


__all__ = ["get_ref_file_args", "run_extract1d"]


@dataclass
class DetectorModelOrder:
    """Model of detector properties for the ATOCA algorithm for a given spectral order."""

    spectral_order: int
    wavemap: np.ndarray
    spectrace: np.ndarray
    specprofile: np.ndarray
    throughput: Callable
    kernel: WebbKernel | None | np.ndarray
    subarray: str | None = None

    @lazyproperty
    def trace(self):
        """
        Get the x, y, wavelength of the trace after applying the transform.

        Returns
        -------
        xtrace, ytrace, wavetrace : array[float]
            The x, y and wavelength of the trace.
        """
        order_idx = self.spectral_order - 1
        spectrace = self.spectrace
        xtrace = np.arange(CUTOFFS[order_idx])

        # CubicSpline requires monotonically increasing x arr
        if spectrace[0][0] - spectrace[0][1] > 0:
            spectrace = np.flip(spectrace, axis=1)

        trace_interp_y = CubicSpline(spectrace[0], spectrace[1])
        trace_interp_wave = CubicSpline(spectrace[0], spectrace[2])
        ytrace = trace_interp_y(xtrace)
        wavetrace = trace_interp_wave(xtrace)
        return xtrace, ytrace, wavetrace

    @lazyproperty
    def native_grid(self):
        """
        Make a 1d-grid of the pixels boundary based on the wavelength solution.

        Returns
        -------
        wave : array[float]
            Grid of the pixels boundaries at the native sampling (1d array)
        col : array[int]
            The column number of the pixel
        """
        # From wavelength solution
        col, _, wave = self.trace

        # Keep only valid solution ...
        idx_valid = np.isfinite(wave)
        # ... and should correspond to subsequent columns
        is_subsequent = np.diff(col[idx_valid]) == 1
        if not is_subsequent.all():
            msg = f"Wavelength solution for order {self.spectral_order} contains gaps."
            log.warning(msg)
        wave = wave[idx_valid]
        col = col[idx_valid]
        log.debug(f"Wavelength range for order {self.spectral_order}: ({wave[[0, -1]]})")

        # Sort
        idx_sort = np.argsort(wave)
        wave = wave[idx_sort]
        col = col[idx_sort]

        return wave, col

    def get_grid_from_trace(self, n_os):
        """
        Make a 1d-grid of the pixels boundary based on the wavelength solution.

        Parameters
        ----------
        n_os : int or array
            The oversampling factor of the wavelength grid used when solving for
            the uncontaminated flux.

        Returns
        -------
        array[float]
            Grid of the pixels boundaries at the native sampling (1d array)
        """
        wave, _ = self.native_grid

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
    detector_models = []
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

        # make throughput interpolator
        thru = throughput_soss(
            pastasoss_ref.throughputs[order_idx].wavelength[:],
            pastasoss_ref.throughputs[order_idx].throughput[:],
        )

        detector_model = DetectorModelOrder(
            spectral_order=order,
            subarray=ref_files["subarray"],
            wavemap=wavemap,
            specprofile=specprofile,
            throughput=thru,
            kernel=None,
            spectrace=spectraces[order_idx],
        )

        # Build a kernel for this order
        wv_cent = np.zeros(wavemap.shape[1])

        # Get central wavelength as a function of columns
        col, row, wv = detector_model.trace
        wv_cent[col] = wv

        # Set invalid values to zero
        idx_invalid = ~np.isfinite(wv_cent)
        wv_cent[idx_invalid] = 0.0

        kernel = WebbKernel(speckernel_ref.wavelengths, speckernel_ref.kernels, wv_cent, n_pix)
        valid_wavemap = (speckernel_wv_range[0] <= wavemap) & (wavemap <= speckernel_wv_range[1])
        wavemap = np.where(valid_wavemap, wavemap, 0.0)
        detector_model.kernel = kernel
        detector_models.append(detector_model)

    return detector_models


def _estim_flux_first_order(
    scidata_bkg, scierr, scimask, order1_model, order2_profile, mask_trace_profile, threshold=1e-4
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
    order1_model : DetectorModelOrder
        The model for order 1.
    order2_profile : array
        The spatial profile for order 2.
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
    # Unpack ref_file arguments
    wave_map = order1_model.wavemap
    spat_pro = order1_model.specprofile
    thrpt = order1_model.throughput

    # Define wavelength grid based on order 1 only (so first index)
    wave_grid = grid_from_map_with_extrapolation(wave_map, spat_pro, n_os=1)

    # Mask parts contaminated by order 2 based on its spatial profile
    mask = (order2_profile >= threshold) | mask_trace_profile[0]

    # Init extraction without convolution kernel (so extract the spectrum at order 1 resolution)
    engine = ExtractionEngine(
        [wave_map], [spat_pro], [thrpt], [None], wave_grid, [mask], global_mask=scimask, orders=[1]
    )

    # Extract estimate
    spec_estimate = engine(scidata_bkg, scierr)

    # Interpolate
    idx = np.isfinite(spec_estimate)
    return UnivariateSpline(wave_grid[idx], spec_estimate[idx], k=3, s=0, ext=0)


def _make_decontamination_grid(detector_models, rtol, max_grid_size, estimate, n_os):
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
    detector_models : list[DetectorModelOrder]
        The list of detector models, one per spectral order.
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
    spectral_orders = [1, 2]
    grids_ord = {}
    for sp_ord in spectral_orders:
        order_model = detector_models[sp_ord - 1]
        grids_ord[sp_ord] = order_model.get_grid_from_trace(n_os)

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


def _f_to_spec(f_order, grid_order, detector_model, pixel_grid, mask):
    """
    Bin the flux to the pixel grid and build a SpecModel.

    Parameters
    ----------
    f_order : np.array
        The solution f_k of the linear system.
    grid_order : np.array
        The wavelength grid of the solution, usually oversampled compared to the pixel grid.
    detector_model : DetectorModelOrder
        The model for the spectral order.
    pixel_grid : np.array
        The pixel grid to which the flux should be binned.
    mask : np.array
        The mask of the pixels to be extracted.

    Returns
    -------
    spec : SpecModel
        The SpecModel containing the extracted spectrum.
    """
    # Build 1d spectrum integrated over pixels
    pixel_grid = pixel_grid[np.newaxis, :]
    model = ExtractionEngine(
        [pixel_grid],
        [np.ones_like(pixel_grid)],
        [detector_model.throughput],
        [detector_model.kernel],
        wave_grid=grid_order,
        mask_trace_profile=[mask[np.newaxis, :]],
        orders=[detector_model.spectral_order],
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
    spec.spectral_order = detector_model.spectral_order

    return spec


def _build_tracemodel_order(engine, detector_model, f_k, mask):
    # Take only the order's specific ref_files
    # And give the identity kernel to the Engine (so no convolution)
    i_order = detector_model.spectral_order - 1

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
    model = ExtractionEngine(
        [detector_model.wavemap],
        [detector_model.specprofile],
        [detector_model.throughput],
        [np.array([1.0])],
        wave_grid=grid_order,
        mask_trace_profile=[mask],
        orders=[sp_ord],
    )

    # Project on detector and save in dictionary
    tracemodel_ord = model.rebuild(flux_order, fill_value=np.nan)

    # Build 1d spectrum integrated over pixels
    pixel_wave_grid, valid_cols = detector_model.native_grid
    spec_ord = _f_to_spec(
        flux_order,
        grid_order,
        detector_model,
        pixel_wave_grid,
        np.all(mask, axis=0)[valid_cols],
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


def _model_image(
    scidata_bkg,
    scierr,
    scimask,
    refmask,
    detector_models,
    box_weights,
    order_list,
    tikfac=None,
    threshold=1e-4,
    n_os=2,
    wave_grid=None,
    estimate=None,
    rtol=1e-3,
    max_grid_size=1000000,
):
    """
    Perform the spectral extraction on a single image.

    The extraction takes place in three stages.
    First, Order 1 is extracted together with the long-wave portion of Order 2,
    since these are overlapping.
    Second, the short-wave portion of Order 2 is extracted standalone, since it
    is assumed not to overlap much with the other orders.
    Third, all of Order 3 (if requested) is extracted standalone, since it too is assumed
    not to overlap with the other orders.

    Parameters
    ----------
    scidata_bkg : array[float]
        A single background subtracted NIRISS SOSS detector image.
    scierr : array[float]
        The uncertainties corresponding to the detector image.
    scimask : array[bool]
        Pixel mask to apply to detector image.
    refmask : array[bool]
        Pixels that should never be reconstructed e.g. the reference pixels.
    detector_models : list[DetectorModelOrder]
        A list of models for each spectral order.
    box_weights : dict
        A dictionary of the weights (for each order) used in the box extraction.
        The weights for each order are 2d arrays with the same size as the detector.
    order_list : list
        List of spectral orders to extract. Must be either [1,2] or [1,2,3].
    tikfac : float, optional
        The Tikhonov regularization factor used when solving for
        the uncontaminated flux. If not specified, the optimal Tikhonov factor
        is calculated.
    threshold : float
        The threshold value for using pixels based on the spectral profile.
        Default value is 1e-4.
    n_os : int, optional
        The oversampling factor of the wavelength grid used when solving for
        the uncontaminated flux. If not specified, defaults to 2.
    wave_grid : ndarray, optional
        Wavelength grid used by ATOCA to model each pixel valid pixel of the detector.
        If not given, the grid is determined based on an estimate of the flux (estimate),
        the relative tolerance (rtol) required on each pixel model and
        the maximum grid size (max_grid_size).
    estimate : UnivariateSpline or None
         Estimate of the target flux as a function of wavelength in microns.
    rtol : float
        The relative tolerance needed on a pixel model. It is used to determine the sampling
        of wave_grid when the input wave_grid is None. Default is 1e-3.
    max_grid_size : int
        Maximum grid size allowed when wave_grid is None.
        Default is 1000000.

    Returns
    -------
    tracemodels : dict
        Dictionary of the modeled detector images for each order.
    tikfac : float
        Optimal Tikhonov factor used in extraction
    logl : float
        Log likelihood value associated with the Tikhonov factor selected.
    wave_grid : 1d array
        The wavelengths at which the spectra were extracted. Same as wave_grid
        if specified as input.
    spec_list : list of SpecModel
        List of the underlying spectra for each integration and order.
        The tikhonov tests are also included.
    """
    order_strs = [f"Order {order}" for order in order_list]
    order_indices = [order - 1 for order in order_list]

    # Some error values are 0, we need to mask those pixels for the extraction engine.
    scimask = scimask | ~(scierr > 0)

    # Define mask based on box aperture
    # (we want to model each contaminated pixels that will be extracted)
    mask_trace_profile = [(~(box_weights[order_strs[i]] > 0)) | (refmask) for i in order_indices]

    # Define mask of pixel to model (all pixels inside box aperture)
    global_mask = np.all(mask_trace_profile, axis=0).astype(bool)

    # Rough estimate of the underlying flux
    if (tikfac is None or wave_grid is None) and estimate is None:
        order_2_profile = detector_models[1].specprofile
        estimate = _estim_flux_first_order(
            scidata_bkg, scierr, scimask, detector_models[0], order_2_profile, mask_trace_profile
        )

    # Generate grid based on estimate if not given
    if wave_grid is None:
        log.info(f"wave_grid not given: generating grid based on rtol={rtol}")
        wave_grid = _make_decontamination_grid(
            detector_models[:2], rtol, max_grid_size, estimate, n_os
        )
        log.debug(
            f"wave_grid covering from {wave_grid.min()} to {wave_grid.max()}"
            f" with {wave_grid.size} points"
        )
    else:
        log.info("Using previously computed or user specified wavelength grid.")

    # Initialize the Engine for combined extraction of orders 1 and 2
    engine = ExtractionEngine(
        [model.wavemap for model in detector_models[:2]],
        [model.specprofile for model in detector_models[:2]],
        [model.throughput for model in detector_models[:2]],
        [model.kernel for model in detector_models[:2]],
        wave_grid=wave_grid,
        mask_trace_profile=mask_trace_profile[:2],
        global_mask=scimask,
        threshold=threshold,
        orders=[1, 2],
    )

    spec_list = []
    if tikfac is None:
        log.info("Solving for the optimal Tikhonov factor.")
        save_tiktests = True

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
        for order in [1, 2]:
            for idx in range(len(all_tests["factors"])):
                f_k = all_tests["solution"][idx, :]
                args = (engine, detector_models[order - 1], f_k, global_mask)
                _, spec_ord = _build_tracemodel_order(*args)
                _populate_tikho_attr(spec_ord, all_tests, idx, order)
                spec_ord.meta.soss_extract1d.color_range = "RED"

                # Add the result to spec_list
                spec_list.append(spec_ord)
    else:
        save_tiktests = False

    log.info(f"Using a Tikhonov factor of {tikfac}")

    # Run the extract method of the Engine.
    f_k = engine(scidata_bkg, scierr, tikhonov=True, factor=tikfac)

    # Compute the log-likelihood of the best fit.
    logl = engine.compute_likelihood(f_k, scidata_bkg, scierr)

    log.info(f"Optimal solution has a log-likelihood of {logl}")

    # Create a new instance of the engine for evaluating the trace model.
    # This allows bad pixels and pixels below the threshold to be reconstructed as well.
    # Model the traces for each order separately.
    tracemodels = {}
    for i_order, order in enumerate([1, 2]):
        log.debug(f"Building the model image of {order}.")

        args = (engine, detector_models[i_order], f_k, global_mask)
        tracemodel_ord, spec_ord = _build_tracemodel_order(*args)
        spec_ord.meta.soss_extract1d.factor = tikfac
        spec_ord.meta.soss_extract1d.color_range = "RED"
        spec_ord.meta.soss_extract1d.type = "OBSERVATION"

        # Project on detector and save in dictionary
        tracemodels[order_strs[i_order]] = tracemodel_ord

        # Add the result to spec_list
        spec_list.append(spec_ord)

    # Make a null tracemodel to be overwritten later for order 3
    if 3 in order_list:
        tracemodels["Order 3"] = np.zeros_like(tracemodels["Order 2"]) * np.nan

    # Model the blue part of order 2 and all of order 3 assuming they are well-separated
    # from order 1
    for order in order_list[1:]:
        order_model = detector_models[order - 1]
        if order_model.subarray == "SUBSTRIP96":
            continue
        if order not in order_list:
            continue
        if order == 2:
            cutoff = ORDER2_SEPARATION_CUTOFF[1]
        else:
            cutoff = None

        idx_order = np.array(order_indices)[np.array(order_list) == order][0]
        order_str = order_strs[idx_order]
        log.info(f"Generate model for well-separated part of {order_str}")

        # Mask for the fit. All valid pixels inside box aperture
        mask_fit = mask_trace_profile[idx_order] | scimask

        # Build 1d spectrum integrated over pixels
        pixel_wave_grid, valid_cols = order_model.native_grid

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
                scidata_bkg,
                scierr,
                order_model,
                mask_fit,
                global_mask,
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
            model = np.nan * np.ones_like(scidata_bkg)

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


def _compute_box_weights(order_models, shape, width, orders_requested):
    """
    Determine the weights for the box extraction.

    Parameters
    ----------
    order_models : list[DetectorModelOrder]
        Models of the detector and trace properties, one per spectral order.
    shape : tuple
        The shape of the detector image.
    width : int
        The width of the box aperture.
    orders_requested : list[int]
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
        order_idx = order_integer - 1

        log.debug(f"Compute box weights for {order}.")

        # Define the box aperture
        xtrace, ytrace, wavelengths[order] = order_models[order_idx].trace
        box_weights[order] = get_box_weights(ytrace, width, shape, cols=xtrace)

    return box_weights, wavelengths


def _decontaminate_image(scidata_bkg, tracemodels, subarray):
    """
    Perform decontamination of the image based on the trace models.

    Parameters
    ----------
    scidata_bkg : array
        A single background subtracted NIRISS SOSS detector image.
    tracemodels : dict
        Dictionary of the modeled detector images for each order.
    subarray : str
        The subarray used for the observation.

    Returns
    -------
    decontaminated_data : dict
        Dictionary of the decontaminated data for each order.
    """
    # Which orders to extract.
    if subarray == "SUBSTRIP96":
        order_list = [1, 2]
    else:
        order_list = [1, 2, 3]

    order_str = {order: f"Order {order}" for order in order_list}

    # List of modeled orders
    mod_order_list = tracemodels.keys()

    # Create dictionaries for the output images.
    decontaminated_data = {}

    log.debug("Performing the decontamination.")

    # Extract each order from order list
    for order_integer in order_list:
        # Order string-name is used more often than integer-name
        order = order_str[order_integer]

        # Decontaminate using all other modeled orders
        decont = scidata_bkg
        for mod_order in mod_order_list:
            if mod_order != order:
                log.debug(f"Decontaminating {order} from {mod_order} using model.")
                is_valid = np.isfinite(tracemodels[mod_order])
                decont = decont - np.where(is_valid, tracemodels[mod_order], 0.0)

        # Save results
        decontaminated_data[order] = decont

    return decontaminated_data


def _model_single_order(
    data_order,
    err_order,
    detector_model,
    mask_fit,
    mask_rebuild,
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
    detector_model : DetectorModelOrder
        The model for the spectral order to be extracted.
    mask_fit : np.array
        Mask determining the aperture used for extraction. This typically includes
        detector bad pixels and any pixels that are not part of the trace
    mask_rebuild : np.array
        Mask determining the aperture used for rebuilding the trace. This typically includes
        only pixels that do not belong to either spectral trace, i.e., regions of the detector
        where no real data could exist.
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
    order = detector_model.spectral_order

    # The throughput and kernel is not needed here
    # set them so they have no effect on the extraction.
    def throughput(wavelength):
        return np.ones_like(wavelength)

    kernel = np.array([1.0])

    # Define wavelength grid with oversampling of 3 (should be enough)
    wave_grid_os = oversample_grid(wave_grid, n_os=3)

    # Initialize the Engine.
    engine = ExtractionEngine(
        [detector_model.wavemap],
        [detector_model.specprofile],
        [throughput],
        [kernel],
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
                detector_model,
                wave_grid,
                np.all(mask_rebuild, axis=0)[valid_cols],
            )
            _populate_tikho_attr(spec_ord, all_tests, idx, order)

            # Add the result to spec_list
            spec_list.append(spec_ord)

    # Rebuild trace, including bad pixels
    engine = ExtractionEngine(
        [detector_model.wavemap],
        [detector_model.specprofile],
        [throughput],
        [kernel],
        wave_grid=wave_grid_os,
        mask_trace_profile=[mask_rebuild],
        orders=[order],
    )
    model = engine.rebuild(f_k_final, fill_value=np.nan)

    # Build 1d spectrum integrated over pixels
    spec_ord = _f_to_spec(
        f_k_final,
        wave_grid_os,
        detector_model,
        wave_grid,
        np.all(mask_rebuild, axis=0)[valid_cols],
    )
    spec_ord.meta.soss_extract1d.factor = tikfac
    spec_ord.meta.soss_extract1d.type = "OBSERVATION"

    # Add the result to spec_list
    spec_list.append(spec_ord)
    return model, spec_list


# Remove bad pixels that are not modeled for pixel number
def _extract_image(
    decontaminated_data, scierr, scimask, box_weights, bad_pix="model", tracemodels=None
):
    """
    Perform the box-extraction on the image using the trace model to correct for contamination.

    Parameters
    ----------
    decontaminated_data : array[float]
        A single background subtracted NIRISS SOSS detector image.
    scierr : array[float]
        The uncertainties corresponding to the detector image.
    scimask : array[float]
        Pixel mask to apply to the detector image.
    box_weights : dict
        A dictionary of the weights (for each order) used in the box extraction.
        The weights for each order are 2d arrays with the same size as the detector.
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
    order_list = box_weights.keys()

    # Create dictionaries for the output spectra.
    fluxes, fluxerrs, npixels = {}, {}, {}

    log.info("Performing the box extraction.")

    # Extract each order from order list
    for order in order_list:
        log.debug(f"Extracting {order}.")

        # Define the box aperture
        box_w_ord = box_weights[order]

        # Decontaminate using all other modeled orders
        decont = decontaminated_data[order]

        # Deal with bad pixels if required.
        if bad_pix == "model":
            # Model the bad pixels decontaminated image when available
            try:
                # Some pixels might not be modeled by the bad pixel models
                is_modeled = np.isfinite(tracemodels[order])
                # Replace bad pixels
                decont = np.where(scimask & is_modeled, tracemodels[order], decont)

                log.debug(f"Bad pixels in {order} are replaced with trace model.")

                # Replace error estimate of the bad pixels
                # using other valid pixels of similar value.
                # The pixel to be estimated are the masked pixels in the region of extraction
                # with available model.
                extraction_region = box_w_ord > 0
                pix_to_estim = extraction_region & scimask & is_modeled
                # Use only valid pixels (not masked) in the extraction region
                # for the empirical estimation
                valid_pix = extraction_region & ~scimask
                scierr_ord = estim_error_nearest_data(scierr, decont, pix_to_estim, valid_pix)

                # Update the scimask for box extraction:
                # the pixels that are modeled are not masked anymore, so set to False.
                # Note that they have to be in the extraction region
                # to ensure that scierr is also valid
                scimask_ord = np.where(is_modeled, False, scimask)

            except KeyError:
                # Keep same mask and error
                scimask_ord = scimask
                scierr_ord = scierr
                log.warning(
                    f"Bad pixels in {order} will be masked instead of modeled: "
                    "trace model unavailable."
                )
        else:
            scimask_ord = scimask
            scierr_ord = scierr
            log.info(f"Bad pixels in {order} will be masked.")

        # Perform the box extraction and save
        out = box_extract(decont, scierr_ord, scimask_ord, box_w_ord)
        _, fluxes[order], fluxerrs[order], npixels[order] = out

    return fluxes, fluxerrs, npixels


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
    order_str_to_int = {f"Order {order}": order for order in order_list}

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

    all_tracemodels, all_box_weights = {}, {}

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
    order_models = get_ref_file_args(ref_files, orders_requested=order_list)
    # Pre-compute the weights for box extraction (used in modeling and extraction)
    shape = (cube_model.data.shape[1], cube_model.data.shape[2])
    box_weights, wavelengths = _compute_box_weights(
        order_models,
        shape,
        width=soss_kwargs["width"],
        orders_requested=order_list,
    )

    # Loop over images.
    output_spec_list = {}
    for i in range(nimages):
        log.info(f"Processing integration {i + 1} of {nimages}.")

        # Unpack the i-th image, set dtype to float64 and convert DQ to boolean mask.
        scidata = cube_model.data[i].astype("float64")
        scierr = cube_model.err[i].astype("float64")
        scimask = np.bitwise_and(cube_model.dq[i], dqflags.pixel["DO_NOT_USE"]).astype(bool)
        refmask = bitfield_to_boolean_mask(
            cube_model.dq[i], ignore_flags=dqflags.pixel["REFERENCE_PIXEL"], flip_bits=True
        )

        # Make sure there aren't any nans not flagged in scimask
        not_finite = ~(np.isfinite(scidata) & np.isfinite(scierr))
        if (not_finite & ~scimask).any():
            log.warning(
                "Input contains invalid values that "
                "are not flagged correctly in the dq map. "
                "They will be masked for the following procedure."
            )
            scimask |= not_finite
            refmask &= ~not_finite

        # Perform background correction.
        if soss_kwargs["subtract_background"]:
            log.info("Applying background subtraction.")
            bkg_mask = make_background_mask(scidata, width=40)
            scidata_bkg, col_bkg = soss_background(scidata, scimask, bkg_mask)
        else:
            log.info("Skip background subtraction.")
            scidata_bkg = scidata
            col_bkg = np.zeros(scidata.shape[1])

        # FIXME: hardcoding the substrip96 weights to unity is a band-aid solution
        if subarray == "SUBSTRIP96":
            box_weights["Order 2"] = np.ones((96, 2048))

        # Model the traces based on optics filter configuration (CLEAR or F277W)
        if soss_filter == "CLEAR" and generate_model:
            # Model the image.
            kwargs = {}
            kwargs["order_list"] = order_list
            kwargs["estimate"] = estimate
            kwargs["tikfac"] = soss_kwargs["tikfac"]
            kwargs["max_grid_size"] = soss_kwargs["max_grid_size"]
            kwargs["rtol"] = soss_kwargs["rtol"]
            kwargs["n_os"] = soss_kwargs["n_os"]
            kwargs["wave_grid"] = wave_grid
            kwargs["threshold"] = soss_kwargs["threshold"]

            result = _model_image(
                scidata_bkg, scierr, scimask, refmask, order_models, box_weights, **kwargs
            )
            tracemodels, soss_kwargs["tikfac"], _, wave_grid, spec_list = result

            # Add atoca spectra to multispec for output
            for spec in spec_list:
                # If it was a test, not the best spectrum,
                # int_num is already set to 0.
                if not hasattr(spec, "int_num"):
                    spec.int_num = i + 1
                output_atoca.spec.append(spec)

        elif soss_filter != "CLEAR" and generate_model:
            # No model can be fit for F277W yet, missing throughput reference files.
            msg = f"No extraction possible for filter {soss_filter}."
            log.critical(msg)
            raise ValueError(msg)
        else:
            # Return empty tracemodels
            tracemodels = {}

        # Decontaminate the data using trace models (if tracemodels not empty)
        data_to_extract = _decontaminate_image(scidata_bkg, tracemodels, subarray)

        if soss_kwargs["bad_pix"] == "model":
            # Generate new trace models for each individual decontaminated orders
            bad_pix_models = tracemodels
        else:
            bad_pix_models = None

        # Use the bad pixel models to perform a de-contaminated extraction.
        kwargs = {}
        kwargs["bad_pix"] = soss_kwargs["bad_pix"]
        kwargs["tracemodels"] = bad_pix_models
        result = _extract_image(data_to_extract, scierr, scimask, box_weights, **kwargs)
        fluxes, fluxerrs, npixels = result

        # Save trace models for output reference
        for order in tracemodels:
            # Initialize a list for first integration
            if i == 0:
                all_tracemodels[order] = []
            # Put NaNs to zero
            model_ord = tracemodels[order]
            model_ord = np.where(np.isfinite(model_ord), model_ord, 0.0)
            # Save as a list (convert to array at the end)
            all_tracemodels[order].append(model_ord)

        # Save box weights for output reference
        for order in box_weights:
            # Initialize a list for first integration
            if i == 0:
                all_box_weights[order] = []
            all_box_weights[order].append(box_weights[order])
        # Copy spectral data for each order into the output model.
        for order in fluxes.keys():
            table_size = len(wavelengths[order])

            out_table = np.zeros(table_size, dtype=datamodels.SpecModel().spec_table.dtype)
            out_table["WAVELENGTH"] = wavelengths[order][:table_size]
            out_table["FLUX"] = fluxes[order][:table_size]
            out_table["FLUX_ERROR"] = fluxerrs[order][:table_size]
            out_table["DQ"] = np.zeros(table_size)
            out_table["BACKGROUND"] = col_bkg[:table_size]
            out_table["NPIXELS"] = npixels[order][:table_size]

            spec = datamodels.SpecModel(spec_table=out_table)

            # Add integration number and spectral order
            spec.spectral_order = order_str_to_int[order]
            spec.int_num = i + 1  # integration number starts at 1, not 0 like python

            if order in output_spec_list:
                output_spec_list[order].append(spec)
            else:
                output_spec_list[order] = [spec]

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
        order_int = order_str_to_int[order]
        setattr(output_references, f"order{order_int}", tracemod_ord)

    for order in all_box_weights:
        # Convert from list to array
        box_w_ord = np.array(all_box_weights[order])
        # Save
        order_int = order_str_to_int[order]
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
