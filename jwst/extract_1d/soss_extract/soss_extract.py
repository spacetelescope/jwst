import logging

import numpy as np
from scipy.interpolate import UnivariateSpline

from ..extract import populate_time_keywords
from ...lib import pipe_utils
from ... import datamodels
from ...datamodels import dqflags, SossWaveGridModel
from astropy.nddata.bitmask import bitfield_to_boolean_mask

from .soss_syscor import make_background_mask, soss_background
from .soss_solver import solve_transform, transform_wavemap, transform_profile, transform_coords
from .atoca import ExtractionEngine
from .atoca_utils import (ThroughputSOSS, WebbKernel, grid_from_map, mask_bad_dispersion_direction,
                          make_combined_adaptive_grid, get_wave_p_or_m, oversample_grid)
from .soss_boxextract import get_box_weights, box_extract, estim_error_nearest_data

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def get_ref_file_args(ref_files, transform):
    """Prepare the reference files for the extraction engine.
    Parameters
    ----------
    ref_files : dict
        A dictionary of the reference file DataModels.
    transform :  array or list
        A 3-element array describing the rotation and translation to apply
        to the reference files in order to match the observation.

    Returns
    -------
    tuple
        The reference file args used with the extraction engine:
        (wavemaps, specprofiles, throughputs, kernels)
    """

    # The wavelength maps for order 1 and 2.
    wavemap_ref = ref_files['wavemap']

    ovs = wavemap_ref.map[0].oversampling
    pad = wavemap_ref.map[0].padding

    wavemap_o1 = transform_wavemap(transform, wavemap_ref.map[0].data, ovs, pad)
    wavemap_o2 = transform_wavemap(transform, wavemap_ref.map[1].data, ovs, pad)

    # Make sure all pixels follow the expected direction of the dispersion
    wavemap_o1, flag_o1 = mask_bad_dispersion_direction(wavemap_o1)
    wavemap_o2, flag_o2 = mask_bad_dispersion_direction(wavemap_o2)

    # Warn if not all pixels were corrected
    msg_warning = 'Some pixels in order {} do not follow the expected dispersion axis'
    if not flag_o1:
        log.warning(msg_warning.format(1))
    if not flag_o2:
        log.warning(msg_warning.format(2))

    # The spectral profiles for order 1 and 2.
    specprofile_ref = ref_files['specprofile']
    ovs = specprofile_ref.profile[0].oversampling
    pad = specprofile_ref.profile[0].padding

    specprofile_o1 = transform_profile(transform, specprofile_ref.profile[0].data, ovs, pad, norm=False)
    specprofile_o2 = transform_profile(transform, specprofile_ref.profile[1].data, ovs, pad, norm=False)

    # The throughput curves for order 1 and 2.
    spectrace_ref = ref_files['spectrace']

    throughput_o1 = ThroughputSOSS(spectrace_ref.trace[0].data['WAVELENGTH'], spectrace_ref.trace[0].data['THROUGHPUT'])
    throughput_o2 = ThroughputSOSS(spectrace_ref.trace[1].data['WAVELENGTH'], spectrace_ref.trace[1].data['THROUGHPUT'])

    # The spectral kernels.
    speckernel_ref = ref_files['speckernel']
    ovs = speckernel_ref.meta.spectral_oversampling
    n_pix = 2 * speckernel_ref.meta.halfwidth + 1

    # Take the centroid of each trace as a grid to project the WebbKernel
    # WebbKer needs a 2d input, so artificially add axis
    wave_maps = [wavemap_o1, wavemap_o2]
    centroid = dict()
    for wv_map, order in zip(wave_maps, [1, 2]):
        # Needs the same number of columns as the detector. Put zeros where not define.
        wv_cent = np.zeros((1, wv_map.shape[1]))
        # Get central wavelength as a function of columns
        col, _, wv = get_trace_1d(ref_files, transform, order)
        wv_cent[:, col] = wv
        # Set invalid values to zero
        idx_invalid = ~np.isfinite(wv_cent)
        wv_cent[idx_invalid] = 0.0
        centroid[order] = wv_cent

    # Get kernels
    kernels_o1 = WebbKernel(speckernel_ref.wavelengths, speckernel_ref.kernels, centroid[1], ovs, n_pix)
    kernels_o2 = WebbKernel(speckernel_ref.wavelengths, speckernel_ref.kernels, centroid[2], ovs, n_pix)

    # Temporary(?) fix to make sure that the kernels can cover the wavelength maps
    speckernel_wv_range = [np.min(speckernel_ref.wavelengths), np.max(speckernel_ref.wavelengths)]
    valid_wavemap = (speckernel_wv_range[0] <= wavemap_o1) & (wavemap_o1 <= speckernel_wv_range[1])
    wavemap_o1 = np.where(valid_wavemap, wavemap_o1, 0.)
    valid_wavemap = (speckernel_wv_range[0] <= wavemap_o2) & (wavemap_o2 <= speckernel_wv_range[1])
    wavemap_o2 = np.where(valid_wavemap, wavemap_o2, 0.)

    return [wavemap_o1, wavemap_o2], [specprofile_o1, specprofile_o2], [throughput_o1, throughput_o2], [kernels_o1, kernels_o2]


def get_trace_1d(ref_files, transform, order, cols=None):
    """Get the x, y, wavelength of the trace after applying the transform.
    Parameters
    ----------
    ref_files : dict
        A dictionary of the reference file DataModels.
    transform : array or list
        A 3-element list or array describing the rotation and translation
        to apply to the reference files in order to match the
        observation.
    order : int
        The spectral order for which to return the trace parameters.
    cols : array[int], optional
        The columns on the detector for which to compute the trace
        parameters. If not given, all columns will be computed.

    Returns
    -------
    xtrace, ytrace, wavetrace : array[float]
        The x, y and wavelength of the trace.
    """

    if cols is None:
        xtrace = np.arange(2048)
    else:
        xtrace = cols

    spectrace_ref = ref_files['spectrace']

    # Read x, y, wavelength for the relevant order.
    xref = spectrace_ref.trace[order - 1].data['X']
    yref = spectrace_ref.trace[order - 1].data['Y']
    waveref = spectrace_ref.trace[order - 1].data['WAVELENGTH']

    # Rotate and shift the positions based on transform.
    angle, xshift, yshift = transform
    xrot, yrot = transform_coords(angle, xshift, yshift, xref, yref)

    # Interpolate y and wavelength to the requested columns.
    sort = np.argsort(xrot)
    ytrace = np.interp(xtrace, xrot[sort], yrot[sort])
    wavetrace = np.interp(xtrace, xrot[sort], waveref[sort])

    return xtrace, ytrace, wavetrace


def estim_flux_first_order(scidata_bkg, scierr, scimask, ref_file_args, mask_trace_profile, threshold=1e-4):
    """
    Parameters
    ----------
    scidata_bkg : array
        A single background subtracted NIRISS SOSS detector image.
    scierr : array
        The uncertainties corresponding to the detector image.
    scimask : array
        Pixel mask to apply to the detector image.
    ref_file_args : tuple
        A tuple of reference file arguments constructed by get_ref_file_args().
    mask_trace_profile: array[bool]
        Mask determining the aperture used for extraction. Set to False where the pixel should be extracted.
    threshold : float, optional:
        The pixels with an aperture[order 2] > `threshold` are considered contaminated
        and will be masked. Default is 1e-4.
    Returns
    -------
    func
        A spline estimator that provides the underlying flux as a function of wavelength
    """

    # Unpack ref_file arguments
    wave_maps, spat_pros, thrpts, _ = ref_file_args

    # Oversampling of 1 to make sure the solution will be stable
    n_os = 1

    # Define wavelength grid based on order 1 only (so first index)
    wave_grid = grid_from_map(wave_maps[0], spat_pros[0], n_os=n_os)

    # Mask parts contaminated by order 2 based on its spatial profile
    mask = ((spat_pros[1] >= threshold) | mask_trace_profile | scimask)

    # Init extraction without convolution kernel (so extract the spectrum at order 1 resolution)
    ref_file_args = [wave_maps[0]], [spat_pros[0]], [thrpts[0]], [np.array([1.])]
    kwargs = {'wave_grid': wave_grid,
              'orders': [1],
              'mask_trace_profile': [mask]}
    engine = ExtractionEngine(*ref_file_args, **kwargs)

    # Extract estimate
    spec_estimate = engine.__call__(data=scidata_bkg, error=scierr)

    # Interpolate
    idx = np.isfinite(spec_estimate)
    estimate_spl = UnivariateSpline(wave_grid[idx], spec_estimate[idx], k=3, s=0, ext=0)

    return estimate_spl


def _mask_wv_map_centroid_outside(wave_maps, ref_files, transform, y_max, orders=(1, 2)):
    """Patch to mask wv_map when centroid outside
    Parameters
    ----------
    wave_maps : array or list
        Wavelength maps
    ref_files : dict
        A dictionary of the reference file DataModels.
    transform : array or list
        A 3-element list or array describing the rotation and translation
        to apply to the reference files in order to match the
        observation.
    y_max : int
        Max value of column to check against centroid location
    orders : tuple[int], optional
        The spectral orders for each wave_map. If not specified, defaults to (1, 2).

    Returns
    -------
    None
        Modifies wave_maps in place to mask out values where centroid is off the detector.
    """
    for wv_map, order in zip(wave_maps, orders):
        # Get centroid wavelength and y position as a function of columns
        _, y_pos, wv = get_trace_1d(ref_files, transform, order)
        # Find min and max wavelengths with centroid inside of detector
        wv = wv[np.isfinite(y_pos)]
        y_pos = y_pos[np.isfinite(y_pos)]
        idx_inside = (0 <= y_pos) & (y_pos <= y_max)
        wv = wv[idx_inside]
        # Set to zeros (mask) values outside
        mask = np.isfinite(wv_map)
        mask[mask] = (np.min(wv) > wv_map[mask]) | (wv_map[mask] > np.max(wv))
        wv_map[mask] = 0.


def get_native_grid_from_trace(ref_files, transform, spectral_order):
    """
    Make a 1d-grid of the pixels boundary and ready for ATOCA ExtractionEngine,
    based on the wavelength solution.
    Parameters
    ----------
    ref_files: dict
        A dictionary of the reference file DataModels.
    transform: array_like
        A 3-elemnt list or array describing the rotation and
        translation to apply to the reference files in order to match the
        observation.
    spectral_order: int
        The spectral order for which to return the trace parameters.
    pix_center: bool
        If True, use pixel center wavelength value to define the grid.
        If False, use pixel boundaries. Default is False.
    Returns
    -------
    Grid of the pixels boundaries at the native sampling (1d array)
    """

    # From wavelenght solution
    col, _, wave = get_trace_1d(ref_files, transform, spectral_order)

    # Keep only valid solution ...
    idx_valid = np.isfinite(wave)
    # ... and should correspond to subsequent columns
    is_subsequent = (np.diff(col[idx_valid]) == 1)
    if not is_subsequent.all():
        msg = f'Wavelength solution for order {spectral_order} contains gaps.'
        log.warning(msg)
    wave = wave[idx_valid]
    col = col[idx_valid]
    log.debug(f'Wavelength range for order {spectral_order}: ({wave[[0, -1]]})')

    # Sort
    idx_sort = np.argsort(wave)
    wave = wave[idx_sort]
    col = col[idx_sort]

    return wave, col


def get_grid_from_trace(ref_files, transform, spectral_order, n_os=1):
    """
    Make a 1d-grid of the pixels boundary and ready for ATOCA ExtractionEngine,
    based on the wavelength solution.
    Parameters
    ----------
    ref_files: dict
        A dictionary of the reference file DataModels.
    transform: array_like
        A 3-elemnt list or array describing the rotation and
        translation to apply to the reference files in order to match the
        observation.
    spectral_order: int
        The spectral order for which to return the trace parameters.
    Returns
    -------
    Grid of the pixels boundaries at the native sampling (1d array)
    """

    wave, _ = get_native_grid_from_trace(ref_files, transform, spectral_order)

    # Use pixel boundaries instead of the center values
    wv_upper_bnd, wv_lower_bnd = get_wave_p_or_m(wave[None, :])
    # `get_wave_p_or_m` returns 2d array, so keep only 1d
    wv_upper_bnd, wv_lower_bnd = wv_upper_bnd[0], wv_lower_bnd[0]
    # Each upper boundary should correspond the the lower boundary
    # of the following pixel, so only need to add the last upper boundary to complete the grid
    wv_upper_bnd, wv_lower_bnd = np.sort(wv_upper_bnd), np.sort(wv_lower_bnd)
    wave_grid = np.append(wv_lower_bnd, wv_upper_bnd[-1])

    # Oversample as needed
    wave_grid = oversample_grid(wave_grid, n_os=n_os)

    return wave_grid


def make_decontamination_grid(ref_files, transform, rtol, max_grid_size, estimate, n_os, wv_range=None):
    ''' Create the grid use for the simultaneous extraction of order 1 and 2.
    The grid is made by:
    1) requiring that it satifsfies the oversampling n_os
    2) trying to reach the specified tolerance for the spectral range shared between order 1 and 2
    3) trying to reach the specified tolerance in the rest of spectral range
    The max_grid_size overrules steps 2) and 3), so the precision may not be reached if
    the grid size needed is too large.
    '''

    # Build native grid for each  orders.
    spectral_orders = [2, 1]
    grids_ord = dict()
    for sp_ord in spectral_orders:
        grids_ord[sp_ord] = get_grid_from_trace(ref_files, transform, sp_ord, n_os=n_os)

    # Build the list of grids given to make_combined_grid.
    # It must be ordered in increasing priority.
    # 1rst priority: shared wavelengths with order 1 and 2.
    # 2nd priority: remaining red part of order 1
    # 3rd priority: remaining blue part of order 2
    # So, split order 2 in 2 parts, the shared wavelenght and the bluemost part
    is_shared = grids_ord[2] >= np.min(grids_ord[1])
    # Make sure order 1 is not more in the blue than order 2
    cond = grids_ord[1] > np.min(grids_ord[2][is_shared])
    grids_ord[1] = grids_ord[1][cond]
    # And make grid list
    all_grids = [grids_ord[2][is_shared], grids_ord[1], grids_ord[2][~is_shared]]

    # Set wavelength range if not given
    if wv_range is None:
        # Cut order 2 at 0.77 (not smaller than that)
        # because there is no contamination there. Can be extracted afterward.
        # In the red, no cut.
        wv_range = [0.77, np.max(grids_ord[1])]

    # Finally, build the list of corresponding estimates.
    # The estimate for the overlapping part is the order 1 estimate.
    # There is no estimate yet for the blue part of order 2, so give a flat spectrum.
    def flat_fct(wv):
        return np.ones_like(wv)

    all_estimates = [estimate, estimate, flat_fct]

    # Generate the combined grid
    kwargs = dict(rtol=rtol, max_total_size=max_grid_size, max_iter=30, grid_range=wv_range)
    combined_grid = make_combined_adaptive_grid(all_grids, all_estimates, **kwargs)

    return combined_grid


def append_tiktests(test_a, test_b):

    out = dict()

    for key in test_a:
        out[key] = np.append(test_a[key], test_b[key], axis=0)

    return out


def populate_tikho_attr(spec, tiktests, idx, sp_ord):

    spec.spectral_order = sp_ord
    spec.meta.soss_extract1d.type = 'TEST'
    spec.meta.soss_extract1d.chi2 = tiktests['chi2'][idx]
    spec.meta.soss_extract1d.chi2_soft_l1 = tiktests['chi2_soft_l1'][idx]
    spec.meta.soss_extract1d.chi2_cauchy = tiktests['chi2_cauchy'][idx]
    spec.meta.soss_extract1d.reg = np.nansum(tiktests['reg'][idx] ** 2)
    spec.meta.soss_extract1d.factor = tiktests['factors'][idx]
    spec.int_num = 0

    return


def f_to_spec(f_order, grid_order, ref_file_args, pixel_grid, mask, sp_ord=0):

    # Make sure the input is not modified
    ref_file_args = ref_file_args.copy()

    # Build 1d spectrum integrated over pixels
    pixel_grid = pixel_grid[np.newaxis, :]
    ref_file_args[0] = [pixel_grid]  # Wavelength map
    ref_file_args[1] = [np.ones_like(pixel_grid)]  # No spatial profile
    model = ExtractionEngine(*ref_file_args,
                             wave_grid=grid_order,
                             mask_trace_profile=[mask[np.newaxis, :]],
                             orders=[sp_ord])
    f_binned = model.rebuild(f_order, fill_value=np.nan)

    pixel_grid = np.squeeze(pixel_grid)
    f_binned = np.squeeze(f_binned)
    # Remove Nans to save space
    is_valid = np.isfinite(f_binned)
    table_size = np.sum(is_valid)
    out_table = np.zeros(table_size, dtype=datamodels.SpecModel().spec_table.dtype)
    out_table['WAVELENGTH'] = pixel_grid[is_valid]
    out_table['FLUX'] = f_binned[is_valid]
    spec = datamodels.SpecModel(spec_table=out_table)
    spec.spectral_order = sp_ord

    return spec


def _build_tracemodel_order(engine, ref_file_args, f_k, i_order, mask, ref_files, transform):

    # Take only the order's specific ref_files
    ref_file_order = [[ref_f[i_order]] for ref_f in ref_file_args]

    # Pre-convolve the extracted flux (f_k) at the order's resolution
    # so that the convolution matrix must not be re-computed.
    flux_order = engine.kernels[i_order].dot(f_k)

    # Then must take the grid after convolution (smaller)
    grid_order = engine.wave_grid_c(i_order)

    # Keep only valid values to make sure there will be no Nans in the order model
    idx_valid = np.isfinite(flux_order)
    grid_order, flux_order = grid_order[idx_valid], flux_order[idx_valid]

    # And give the identity kernel to the Engine (so no convolution)
    ref_file_order[3] = [np.array([1.])]

    # Spectral order
    sp_ord = i_order + 1

    # Build model of the order
    model = ExtractionEngine(*ref_file_order,
                             wave_grid=grid_order,
                             mask_trace_profile=[mask],
                             orders=[sp_ord])

    # Project on detector and save in dictionary
    tracemodel_ord = model.rebuild(flux_order, fill_value=np.nan)

    # Build 1d spectrum integrated over pixels
    pixel_wave_grid, valid_cols = get_native_grid_from_trace(ref_files, transform, sp_ord)
    spec_ord = f_to_spec(flux_order, grid_order, ref_file_order, pixel_wave_grid,
                         np.all(mask, axis=0)[valid_cols], sp_ord=sp_ord)

    return tracemodel_ord, spec_ord


def model_image(scidata_bkg, scierr, scimask, refmask, ref_files, box_weights, subarray, transform=None,
                tikfac=None, threshold=1e-4, n_os=2, wave_grid=None,
                estimate=None, rtol=1e-3, max_grid_size=1000000):
    """Perform the spectral extraction on a single image.

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
    ref_files : dict
        A dictionary of the reference file DataModels.
    box_weights : dict
        A dictionary of the weights (for each order) used in the box extraction.
        The weights for each order are 2d arrays with the same size as the detector.
    subarray : str
        Subarray on which the data were recorded; one of 'SUBSTRIPT96',
        'SUBSTRIP256' or 'FULL'.
    transform : array or list, optional
        A 3-element list or array describing the rotation and translation to
        apply to the reference files in order to match the observation. If not
        specified, the transformation is computed.
    tikfac : float, optional
        The Tikhonov regularization factor used when solving for
        the uncontaminated flux. If not specified, the optimal Tikhonov factor
        is calculated.
    n_os : int, optional
        The oversampling factor of the wavelength grid used when solving for
        the uncontaminated flux. If not specified, defaults to 5.
    threshold : float
        The threshold value for using pixels based on the spectral profile.
        Default value is 1e-4.
    wave_grid : str or SossWaveGridModel or None
        Filename of reference file or SossWaveGridModel containing the wavelength grid used by ATOCA
        to model each pixel valid pixel of the detector. If not given, the grid is determined
        based on an estimate of the flux (estimate), the relative tolerance (rtol)
        required on each pixel model and the maximum grid size (max_grid_size).
    estimate : UnivariateSpline or None
         Estimate of the target flux as a function of wavelength in microns.
    rtol : float
        The relative tolerance needed on a pixel model. It is used to determine the sampling
        of the soss_wave_grid when not directly given. Default is 1e-3.
    max_grid_size : int
        Maximum grid size allowed. It is used when soss_wave_grid is not directly
        to make sure the computation time or the memory used stays reasonable.
        Default is 1000000

    Returns
    -------
    tracemodels : dict
        Dictionary of the modeled detector images for each order.
    tikfac : float
        Optimal Tikhonov factor used in extraction
    logl : float
        Log likelihood value associated with the Tikhonov factor selected.
    wave_grid : 1d array
        Same as wave_grid input
    spec_list : list of SpecModel
        List of the underlying spectra for each integration and order.
        The tikhonov tests are also included.
    """

    # Init list of atoca 1d spectra
    spec_list = []

    # Orders to simulate
    order_list = ['Order 1', 'Order 2']

    # Prepare the reference file arguments.
    ref_file_args = get_ref_file_args(ref_files, transform)

    # Some error values are 0, we need to mask those pixels for the extraction engine.
    scimask = scimask | ~(scierr > 0)

    # Define mask based on box aperture (we want to model each contaminated pixels that will be extracted)
    mask_trace_profile = [~(box_weights[order] > 0) for order in order_list]

    # Define mask of pixel to model (all pixels inside box aperture)
    global_mask = np.all(mask_trace_profile, axis=0) | refmask

    # Rough estimate of the underlying flux
    # Note: estim_flux func is not strictly necessary and factors could be a simple logspace -
    #       dq mask caused issues here and this may need a try/except wrap.
    #       Dev suggested np.logspace(-19, -10, 10)
    if (tikfac is None or wave_grid is None) and estimate is None:
        estimate = estim_flux_first_order(scidata_bkg, scierr, scimask,
                                          ref_file_args, mask_trace_profile[0])

    # Generate grid based on estimate if not given
    if wave_grid is None:
        log.info(f'wave_grid not given: generating grid based on rtol={rtol}')
        wave_grid = make_decontamination_grid(ref_files, transform, rtol, max_grid_size, estimate, n_os)
        log.debug(f'wave_grid covering from {wave_grid.min()} to {wave_grid.max()}')
    else:
        log.info('Using previously computed or user specified wavelength grid.')

#     # Use estimate to evaluate the contribution from each orders to pixels
#     # (Used to determine which pixel to model later)
#     ref_args_estimate = [ref_arg for ref_arg in ref_file_args]
#     # No convolution needed (so give equivalent of identity)
#     ref_args_estimate[3] = [np.array([1.]) for _ in order_list]
#     engine_for_estimate = ExtractionEngine(*ref_args_estimate, wave_grid=wave_grid, mask_trace_profile=mask_trace_profile)
#     models = {order: engine_for_estimate.rebuild(estimate, i_orders=[idx_ord], fill_value=np.nan)
#               for idx_ord, order in enumerate(order_list)}
#     total = np.nansum([models[order] for order in order_list], axis=0)
#     total = np.where((total != 0), total, np.nan)
#     contribution = {order: models[order] / total for order in order_list}

    log.debug('Extracting using transformation parameters {}'.format(transform))

    # Set the c_kwargs using the minimum value of the kernels
    c_kwargs = [{'thresh': webb_ker.min_value} for webb_ker in ref_file_args[3]]

    # Initialize the Engine.
    engine = ExtractionEngine(*ref_file_args,
                              wave_grid=wave_grid,
                              mask_trace_profile=mask_trace_profile,
                              global_mask=scimask,
                              threshold=threshold,
                              c_kwargs=c_kwargs)

    if tikfac is None:

        log.info('Solving for the optimal Tikhonov factor.')

        # Find the tikhonov factor.
        # Initial pass 8 orders of magnitude with 10 grid points.
        guess_factor = engine.estimate_tikho_factors(estimate)
        log_guess = np.log10(guess_factor)
        factors = np.logspace(log_guess - 4, log_guess + 4, 10)
        all_tests = engine.get_tikho_tests(factors, data=scidata_bkg, error=scierr)
        tikfac, mode, _ = engine.best_tikho_factor(tests=all_tests, fit_mode='all')

        # Refine across 4 orders of magnitude.
        tikfac = np.log10(tikfac)
        factors = np.logspace(tikfac - 2, tikfac + 2, 20)
        tiktests = engine.get_tikho_tests(factors, data=scidata_bkg, error=scierr)
        tikfac, mode, _ = engine.best_tikho_factor(tests=tiktests, fit_mode='d_chi2')
        # Add all theses tests to previous ones
        all_tests = append_tiktests(all_tests, tiktests)

        # Save spectra in a list of SingleSpecModels for optional output
        save_tiktests = True
        for i_order, order in enumerate(order_list):
            for idx in range(len(all_tests['factors'])):
                f_k = all_tests['solution'][idx, :]
                args = (engine, ref_file_args, f_k, i_order, global_mask, ref_files, transform)
                _, spec_ord = _build_tracemodel_order(*args)
                populate_tikho_attr(spec_ord, all_tests, idx, i_order + 1)
                spec_ord.meta.soss_extract1d.color_range = 'RED'

                # Add the result to spec_list
                spec_list.append(spec_ord)
    else:
        save_tiktests = False

    log.info('Using a Tikhonov factor of {}'.format(tikfac))

    # Run the extract method of the Engine.
    f_k = engine.__call__(data=scidata_bkg, error=scierr, tikhonov=True, factor=tikfac)

    # Compute the log-likelihood of the best fit.
    logl = engine.compute_likelihood(f_k, same=False)

    log.info('Optimal solution has a log-likelihood of {}'.format(logl))

    # Create a new instance of the engine for evaluating the trace model.
    # This allows bad pixels and pixels below the threshold to be reconstructed as well.
    # Model the order 1 and order 2 trace separately.
    tracemodels = dict()

    for i_order, order in enumerate(order_list):

        log.debug('Building the model image of {}.'.format(order))

        args = (engine, ref_file_args, f_k, i_order, global_mask, ref_files, transform)
        tracemodel_ord, spec_ord = _build_tracemodel_order(*args)
        spec_ord.meta.soss_extract1d.factor = tikfac
        spec_ord.meta.soss_extract1d.color_range = 'RED'
        spec_ord.meta.soss_extract1d.type = 'OBSERVATION'

        # Project on detector and save in dictionary
        tracemodels[order] = tracemodel_ord

        # Add the result to spec_list
        spec_list.append(spec_ord)

    # ###############################
    # Model remaining part of order 2
    # ###############################
    if subarray != 'SUBSTRIP96':
        idx_order2 = 1
        order = idx_order2 + 1
        order_str = 'Order 2'
        log.info('Generate model for blue-most part of order 2')

        # Take only the second order's specific ref_files
        ref_file_order = [[ref_f[idx_order2]] for ref_f in ref_file_args]

        # Mask for the fit. All valid pixels inside box aperture
        mask_fit = mask_trace_profile[idx_order2] | scimask
#         # and extract only what was not already modeled
#         already_modeled = np.isfinite(tracemodels[order_str])
#         mask_fit |= already_modeled

        # Build 1d spectrum integrated over pixels
        pixel_wave_grid, valid_cols = get_native_grid_from_trace(ref_files, transform, order)

        # Hardcode wavelength highest boundary as well.
        # Must overlap with lower limit in make_decontamination_grid
        is_in_wv_range = (pixel_wave_grid < 0.95)
        pixel_wave_grid, valid_cols = pixel_wave_grid[is_in_wv_range], valid_cols[is_in_wv_range]

        # NOTE: This code is currently unused.
        # Remove order 1
        # scidata_order2_decont = scidata_bkg - tracemodels['Order 1']

        # Range of initial tikhonov factors
        tikfac_log_range = np.log10(tikfac) + np.array([-2, 8])

        # Model the remaining part of order 2 with atoca
        model, spec_ord = model_single_order(scidata_bkg, scierr, ref_file_order,
                                             mask_fit, global_mask, order,
                                             pixel_wave_grid, valid_cols, save_tiktests,
                                             tikfac_log_range=tikfac_log_range)

        # Keep only pixels from which order 2 contribution
        # is not already modeled.
        already_modeled = np.isfinite(tracemodels[order_str])
        model = np.where(already_modeled, 0., model)

        # Add to tracemodels
        tracemodels[order_str] = np.nansum([tracemodels[order_str], model], axis=0)

        # Add the result to spec_list
        for sp in spec_ord:
            sp.meta.soss_extract1d.color_range = 'BLUE'
        spec_list += spec_ord

    return tracemodels, tikfac, logl, wave_grid, spec_list


def compute_box_weights(ref_files, transform, subarray, shape, width=40.):

    # Which orders to compute (for modeling, different than extraction).
    if subarray == 'SUBSTRIP96':
        order_list = [1, 2]
    else:
        order_list = [1, 2, 3]

    # Extract each order from order list
    box_weights = dict()
    wavelengths = dict()
    order_str = {order: f'Order {order}' for order in order_list}
    for order_integer in order_list:
        # Order string-name is used more often than integer-name
        order = order_str[order_integer]

        log.debug(f'Compute box weights for order {order}.')

        # Define the box aperture
        xtrace, ytrace, wavelengths[order] = get_trace_1d(ref_files, transform, order_integer)
        box_weights[order] = get_box_weights(ytrace, width, shape, cols=xtrace)

    return box_weights, wavelengths


def decontaminate_image(scidata_bkg, tracemodels, subarray):
    """Perform decontamination of the image based on the trace models"""

    # Which orders to extract.
    if subarray == 'SUBSTRIP96':
        order_list = [1, 2]
    else:
        order_list = [1, 2, 3]

    order_str = {order: f'Order {order}' for order in order_list}

    # List of modeled orders
    mod_order_list = tracemodels.keys()

    # Create dictionaries for the output images.
    decontaminated_data = dict()

    log.debug('Performing the decontamination.')

    # Extract each order from order list
    for order_integer in order_list:
        # Order string-name is used more often than integer-name
        order = order_str[order_integer]

        # Decontaminate using all other modeled orders
        decont = scidata_bkg
        for mod_order in mod_order_list:
            if mod_order != order:
                log.debug(f'Decontaminating {order} from {mod_order} using model.')
                is_valid = np.isfinite(tracemodels[mod_order])
                decont = decont - np.where(is_valid, tracemodels[mod_order], 0.)

        # Save results
        decontaminated_data[order] = decont

    return decontaminated_data


# TODO Add docstring
# TODO Add threshold like in model_image? TO use with the rough (but stable) estimate
def model_single_order(data_order, err_order, ref_file_args, mask_fit,
                       mask_rebuild, order, wave_grid, valid_cols, save_tiktests=False, tikfac_log_range=None):

    # The throughput and kernel is not needed here; set them so they have no effect on the extraction.
    def throughput(wavelength):
        return np.ones_like(wavelength)
    kernel = np.array([1.])

    # Set reference file arguments
    ref_file_args[2] = [throughput]
    ref_file_args[3] = [kernel]

    # ###########################
    # First, generate an estimate
    # (only if the initial guess of tikhonov factor range is not given)
    # ###########################

    if tikfac_log_range is None:
        # Initialize the engine
        engine = ExtractionEngine(*ref_file_args,
                                  wave_grid=wave_grid,
                                  orders=[order],
                                  mask_trace_profile=[mask_fit])

        # Extract estimate
        spec_estimate = engine.__call__(data=data_order, error=err_order)

        # Interpolate
        idx = np.isfinite(spec_estimate)
        estimate_spl = UnivariateSpline(wave_grid[idx], spec_estimate[idx], k=3, s=0, ext=0)

    # ##################################################
    # Second, do the extraction to get the best estimate
    # ##################################################
    # Define wavelength grid with oversampling of 3 (should be enough)
    wave_grid_os = oversample_grid(wave_grid, n_os=3)
    wave_grid_os = wave_grid_os[wave_grid_os > 0.58]

    # Initialize the Engine.
    engine = ExtractionEngine(*ref_file_args,
                              wave_grid=wave_grid_os,
                              orders=[order],
                              mask_trace_profile=[mask_fit])

    # Find the tikhonov factor.
    # Initial pass with tikfac_range.
    if tikfac_log_range is None:
        guess_factor = engine.estimate_tikho_factors(estimate_spl)
        log_guess = np.log10(guess_factor)
        factors = np.log_range(log_guess - 2, log_guess + 8, 10)
    else:
        factors = np.logspace(tikfac_log_range[0], tikfac_log_range[-1] + 8, 10)
    all_tests = engine.get_tikho_tests(factors, data=data_order, error=err_order)
    tikfac, mode, _ = engine.best_tikho_factor(tests=all_tests, fit_mode='all')

    # Refine across 4 orders of magnitude.
    tikfac = np.log10(tikfac)
    factors = np.logspace(tikfac - 2, tikfac + 2, 20)
    tiktests = engine.get_tikho_tests(factors, data=data_order, error=err_order)
    tikfac, mode, _ = engine.best_tikho_factor(tests=tiktests, fit_mode='d_chi2')
    all_tests = append_tiktests(all_tests, tiktests)

    # Run the extract method of the Engine.
    f_k_final = engine.__call__(data=data_order, error=err_order, tikhonov=True, factor=tikfac)

    # Save binned spectra in a list of SingleSpecModels for optional output
    spec_list = []
    if save_tiktests:
        for idx in range(len(all_tests['factors'])):
            f_k = all_tests['solution'][idx, :]

            # Build 1d spectrum integrated over pixels
            spec_ord = f_to_spec(f_k, wave_grid_os, ref_file_args, wave_grid,
                                 np.all(mask_rebuild, axis=0)[valid_cols], sp_ord=order)
            populate_tikho_attr(spec_ord, all_tests, idx, order)

            # Add the result to spec_list
            spec_list.append(spec_ord)

    # ##########################################
    # Third, rebuild trace, including bad pixels
    # ##########################################
    # Initialize the Engine.
    engine = ExtractionEngine(*ref_file_args,
                              wave_grid=wave_grid_os,
                              orders=[order],
                              mask_trace_profile=[mask_rebuild])

    # Project on detector and save in dictionary
    model = engine.rebuild(f_k_final, fill_value=np.nan)

    # Build 1d spectrum integrated over pixels
    spec_ord = f_to_spec(f_k_final, wave_grid_os, ref_file_args, wave_grid,
                         np.all(mask_rebuild, axis=0)[valid_cols], sp_ord=order)
    spec_ord.meta.soss_extract1d.factor = tikfac
    spec_ord.meta.soss_extract1d.type = 'OBSERVATION'

    # Add the result to spec_list
    spec_list.append(spec_ord)

    return model, spec_list


# Remove bad pixels that are not modeled for pixel number
# TODO Update docstring
def extract_image(decontaminated_data, scierr, scimask, box_weights, bad_pix='model', tracemodels=None):
    """Perform the box-extraction on the image, while using the trace model to
    correct for contamination.
    Parameters
    ----------
    decontaminated_data : array[float]
        A single backround subtracted NIRISS SOSS detector image.
    scierr : array[float]
        The uncertainties corresponding to the detector image.
    scimask : array[float]
        Pixel mask to apply to the detector image.
    ref_files : dict
        A dictionary of the reference file DataModels.
    transform : array_like
        A 3-element list or array describing the rotation and translation to
        apply to the reference files in order to match the observation.
    subarray : str
        Subarray on which the data were recorded; one of 'SUBSTRIPT96',
        'SUBSTRIP256' or 'FULL'.
    width : float
        The width of the aperture used to extract the uncontaminated spectrum.
    bad_pix : str
        How to handle the bad pixels. Options are 'masking' and 'model'.
        'masking' will simply mask the bad pixels, such that the number of pixels
        in each column in the box extraction will not be constant, while the
        'model' option uses `tracemodels` to replace the bad pixels.
    tracemodels : dict
        Dictionary of the modeled detector images for each order.
    Returns
    -------
    wavelengths, fluxes, fluxerrs, npixels, box_weights : dict
        Each output is a dictionary, with each extracted order as a key.
    """
    # Init models with an empty dictionary if not given
    if tracemodels is None:
        tracemodels = dict()

    # Which orders to extract (extract the ones with given box aperture).
    order_list = box_weights.keys()

    # Create dictionaries for the output spectra.
    fluxes = dict()
    fluxerrs = dict()
    npixels = dict()

    log.info('Performing the box extraction.')

    # Extract each order from order list
    for order in order_list:

        log.debug(f'Extracting {order}.')

        # Define the box aperture
        box_w_ord = box_weights[order]

        # Decontaminate using all other modeled orders
        decont = decontaminated_data[order]

        # Deal with bad pixels if required.
        if bad_pix == 'model':
            # Model the bad pixels decontaminated image when available
            try:
                # Some pixels might not be modeled by the bad pixel models
                is_modeled = np.isfinite(tracemodels[order])
                # Replace bad pixels
                decont = np.where(scimask & is_modeled, tracemodels[order], decont)

                log.debug(f'Bad pixels in {order} are replaced with trace model.')

                # Replace error estimate of the bad pixels using other valid pixels of similar value.
                # The pixel to be estimated are the masked pixels in the region of extraction
                # with available model.
                extraction_region = (box_w_ord > 0)
                pix_to_estim = (extraction_region & scimask & is_modeled)
                # Use only valid pixels (not masked) in the extraction region for the empirical estimation
                valid_pix = (extraction_region & ~scimask)
                scierr_ord = estim_error_nearest_data(scierr, decont, pix_to_estim, valid_pix)

                # Update the scimask for box extraction:
                # the pixels that are modeled are not masked anymore, so set to False.
                # Note that they have to be in the extraction region to ensure that scierr is also valid
                scimask_ord = np.where(is_modeled, False, scimask)

            except KeyError:
                # Keep same mask and error
                scimask_ord = scimask
                scierr_ord = scierr
                log.warning(f'Bad pixels in {order} will be masked instead of modeled: trace model unavailable.')
        else:
            scimask_ord = scimask
            scierr_ord = scierr
            log.info(f'Bad pixels in {order} will be masked.')

        # Perform the box extraction and save
        out = box_extract(decont, scierr_ord, scimask_ord, box_w_ord)
        _, fluxes[order], fluxerrs[order], npixels[order] = out

    return fluxes, fluxerrs, npixels


def run_extract1d(input_model, spectrace_ref_name, wavemap_ref_name,
                  specprofile_ref_name, speckernel_ref_name, subarray,
                  soss_filter, soss_kwargs):
    """Run the spectral extraction on NIRISS SOSS data.
    Parameters
    ----------
    input_model : DataModel
        The input DataModel.
    spectrace_ref_name : str
        Name of the spectrace reference file.
    wavemap_ref_name : str
        Name of the wavemap reference file.
    specprofile_ref_name : str
        Name of the specprofile reference file.
    speckernel_ref_name : str
        Name of the speckernel reference file.
    subarray : str
        Subarray on which the data were recorded; one of 'SUBSTRIPT96',
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
    generate_model = soss_kwargs['atoca'] or (soss_kwargs['bad_pix'] == 'model')

    # Map the order integer names to the string names
    order_str_2_int = {f'Order {order}': order for order in [1, 2, 3]}

    # Read the reference files.
    spectrace_ref = datamodels.SpecTraceModel(spectrace_ref_name)
    wavemap_ref = datamodels.WaveMapModel(wavemap_ref_name)
    specprofile_ref = datamodels.SpecProfileModel(specprofile_ref_name)
    speckernel_ref = datamodels.SpecKernelModel(speckernel_ref_name)

    ref_files = dict()
    ref_files['spectrace'] = spectrace_ref
    ref_files['wavemap'] = wavemap_ref
    ref_files['specprofile'] = specprofile_ref
    ref_files['speckernel'] = speckernel_ref

    # Initialize the theta, dx, dy transform parameters
    transform = soss_kwargs.pop('transform')
    if transform is None:
        transform = [None, None, None]
    # Save names for logging
    param_name = np.array(['theta', 'x-offset', 'y-offset'])

    # Unpack wave_grid if wave_grid_in was specified.
    wave_grid_in = soss_kwargs['wave_grid_in']
    if wave_grid_in is not None:
        log.info(f'Loading wavelength grid from {wave_grid_in}.')
        wave_grid = datamodels.SossWaveGridModel(wave_grid_in).wavegrid
        # Make sure it as the correct precision
        wave_grid = wave_grid.astype('float64')
    else:
        # wave_grid will be estimated later in the first call of `model_image`
        log.info('Wavelength grid was not specified. Setting `wave_grid` to None.')
        wave_grid = None

    # TODO: Maybe not unpack yet. Use SpecModel attributes
    #       to allow for multiple orders? Create unpacking function.
    # Convert estimate to cubic spline if given.
    # It should be a SpecModel or a file name (string)
    estimate = soss_kwargs.pop('estimate')
    if estimate is not None:
        log.info('Converting the estimate of the flux to spline function.')

        # Convert estimate to cubic spline
        estimate = datamodels.open(estimate)
        wv_estimate = estimate.spec_table['WAVELENGTH']
        flux_estimate = estimate.spec_table['FLUX']
        # Keep only finite values
        idx = np.isfinite(flux_estimate)
        estimate = UnivariateSpline(wv_estimate[idx], flux_estimate[idx], k=3, s=0, ext=0)

    # Initialize the output model.
    output_model = datamodels.MultiSpecModel()
    output_model.update(input_model)  # Copy meta data from input to output.

    # Initialize output spectra returned by ATOCA
    output_atoca = datamodels.MultiSpecModel()
    output_model.update(input_model)

    # Initialize output references (model of the detector and box aperture weights).
    output_references = datamodels.SossExtractModel()
    output_references.update(input_model)

    all_tracemodels = dict()
    all_box_weights = dict()

    # Convert to Cube if datamodels is an ImageModel
    if isinstance(input_model, datamodels.ImageModel):

        cube_model = datamodels.CubeModel(shape=(1, *input_model.shape))
        cube_model.data = input_model.data[None, :, :]
        cube_model.err = input_model.err[None, :, :]
        cube_model.dq = input_model.dq[None, :, :]
        nimages = 1
        log.info('Input is an ImageModel, processing a single integration.')

    elif isinstance(input_model, datamodels.CubeModel):

        cube_model = input_model
        nimages = len(cube_model.data)
        log.info('Input is a CubeModel containing {} integrations.'.format(nimages))

    else:
        msg = "Only ImageModel and CubeModel are implemented for the NIRISS SOSS extraction."
        log.critical(msg)
        raise ValueError(msg)

    # Loop over images.
    for i in range(nimages):

        log.info('Processing integration {} of {}.'.format(i + 1, nimages))

        # Unpack the i-th image, set dtype to float64 and convert DQ to boolean mask.
        scidata = cube_model.data[i].astype('float64')
        scierr = cube_model.err[i].astype('float64')
        scimask = np.bitwise_and(cube_model.dq[i], dqflags.pixel['DO_NOT_USE']).astype(bool)
        refmask = bitfield_to_boolean_mask(cube_model.dq[i], ignore_flags=dqflags.pixel['REFERENCE_PIXEL'],
                                           flip_bits=True)

        # Make sure there aren't any nans not flagged in scimask
        not_finite = ~(np.isfinite(scidata) & np.isfinite(scierr))
        if (not_finite & ~scimask).any():
            log.warning('Input contains invalid values that '
                        'are not flagged correctly in the dq map. '
                        'They will be masked for the following procedure.')
            scimask |= not_finite
            refmask &= ~not_finite

        # Perform background correction.
        if soss_kwargs['subtract_background']:
            log.info('Applying background subtraction.')
            bkg_mask = make_background_mask(scidata, width=40)
            scidata_bkg, col_bkg, npix_bkg = soss_background(scidata, scimask, bkg_mask=bkg_mask)
        else:
            log.info('Skip background subtraction.')
            scidata_bkg = scidata
            col_bkg = np.zeros(scidata.shape[1])

        # Determine the theta, dx, dy transform needed to match scidata trace position to ref file position.
        if None in transform:
            log.info('Solving for the transformation parameters.')

            # Unpack the expected order 1 & 2 positions.
            spectrace_ref = ref_files['spectrace']
            xref_o1 = spectrace_ref.trace[0].data['X']
            yref_o1 = spectrace_ref.trace[0].data['Y']
            xref_o2 = spectrace_ref.trace[1].data['X']
            yref_o2 = spectrace_ref.trace[1].data['Y']

            # Define which parameters to fit
            is_fitted = np.array([value is None for value in transform])

            # Show which parameters are fitted in log
            log.info('Parameters used for fit: ' + ', '.join(param_name[is_fitted]))
            log.info('Fixed parameters: ' + ', '.join(param_name[~is_fitted]))

            # Use the solver on the background subtracted image.
            if subarray == 'SUBSTRIP96' or soss_filter == 'F277W':
                # Use only order 1 to solve theta, dx, dy
                transform = solve_transform(scidata_bkg, scimask, xref_o1, yref_o1,
                                            soss_filter=soss_filter, is_fitted=is_fitted,
                                            guess_transform=transform)
            else:
                transform = solve_transform(scidata_bkg, scimask, xref_o1, yref_o1,
                                            xref_o2, yref_o2, is_fitted=is_fitted,
                                            soss_filter=soss_filter, guess_transform=transform)

        string_list = [f'{name}={value}' for name, value in zip(param_name, transform)]
        log.info('Measured to Reference trace position transform: ' + ', '.join(string_list))

        # Pre-compute the weights for box extraction (used in modeling and extraction)
        args = (ref_files, transform, subarray, scidata_bkg.shape)
        box_weights, wavelengths = compute_box_weights(*args, width=soss_kwargs['width'])

        # Model the traces based on optics filter configuration (CLEAR or F277W)
        if soss_filter == 'CLEAR' and generate_model:

            # Model the image.
            kwargs = dict()
            kwargs['transform'] = transform
            kwargs['estimate'] = estimate
            kwargs['tikfac'] = soss_kwargs['tikfac']
            kwargs['max_grid_size'] = soss_kwargs['max_grid_size']
            kwargs['rtol'] = soss_kwargs['rtol']
            kwargs['n_os'] = soss_kwargs['n_os']
            kwargs['wave_grid'] = wave_grid
            kwargs['threshold'] = soss_kwargs['threshold']

            result = model_image(scidata_bkg, scierr, scimask, refmask, ref_files, box_weights, subarray, **kwargs)
            tracemodels, soss_kwargs['tikfac'], logl, wave_grid, spec_list = result

            # Add atoca spectra to multispec for output
            for spec in spec_list:
                # If it was a test, not the best spectrum,
                # int_num is already set to 0.
                if not hasattr(spec, 'int_num'):
                    spec.int_num = i + 1
                output_atoca.spec.append(spec)

        elif soss_filter != 'CLEAR' and generate_model:
            # No model can be fit for F277W yet, missing throughput reference files.
            msg = f"No extraction possible for filter {soss_filter}."
            log.critical(msg)
            raise ValueError(msg)
        else:
            # Return empty tracemodels and no spec_list
            tracemodels = dict()
            spec_list = None

        # Decontaminate the data using trace models (if tracemodels not empty)
        data_to_extract = decontaminate_image(scidata_bkg, tracemodels, subarray)

        if soss_kwargs['bad_pix'] == 'model':
            # Generate new trace models for each individual decontaminated orders
            # TODO: Use the sum of tracemodels so it can be applied even w/o decontamination
            bad_pix_models = tracemodels
        else:
            bad_pix_models = None

        # Use the bad pixel models to perform a de-contaminated extraction.
        kwargs = dict()
        kwargs['bad_pix'] = soss_kwargs['bad_pix']
        kwargs['tracemodels'] = bad_pix_models
        result = extract_image(data_to_extract, scierr, scimask, box_weights, **kwargs)
        fluxes, fluxerrs, npixels = result

        # Save trace models for output reference
        for order in tracemodels:
            # Initialize a list for first integration
            if i == 0:
                all_tracemodels[order] = []
            # Put NaNs to zero
            model_ord = tracemodels[order]
            model_ord = np.where(np.isfinite(model_ord), model_ord, 0.)
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
            out_table['WAVELENGTH'] = wavelengths[order]
            out_table['FLUX'] = fluxes[order]
            out_table['FLUX_ERROR'] = fluxerrs[order]
            out_table['DQ'] = np.zeros(table_size)
            out_table['BACKGROUND'] = col_bkg
            out_table['NPIXELS'] = npixels[order]

            spec = datamodels.SpecModel(spec_table=out_table)

            # Add integration number and spectral order
            spec.spectral_order = order_str_2_int[order]
            spec.int_num = i + 1  # integration number starts at 1, not 0 like python

            output_model.spec.append(spec)

        output_model.meta.soss_extract1d.width = soss_kwargs['width']
        output_model.meta.soss_extract1d.apply_decontamination = soss_kwargs['atoca']
        output_model.meta.soss_extract1d.tikhonov_factor = soss_kwargs['tikfac']
        output_model.meta.soss_extract1d.delta_x = transform[1]
        output_model.meta.soss_extract1d.delta_y = transform[2]
        output_model.meta.soss_extract1d.theta = transform[0]
        output_model.meta.soss_extract1d.oversampling = soss_kwargs['n_os']
        output_model.meta.soss_extract1d.threshold = soss_kwargs['threshold']
        output_model.meta.soss_extract1d.bad_pix = soss_kwargs['bad_pix']

    # Save output references
    for order in all_tracemodels:
        # Convert from list to array
        tracemod_ord = np.array(all_tracemodels[order])
        # Save
        order_int = order_str_2_int[order]
        setattr(output_references, f'order{order_int}', tracemod_ord)

    for order in all_box_weights:
        # Convert from list to array
        box_w_ord = np.array(all_box_weights[order])
        # Save
        order_int = order_str_2_int[order]
        setattr(output_references, f'aperture{order_int}', box_w_ord)

    if pipe_utils.is_tso(input_model):
        log.info("Populating INT_TIMES keywords from input table.")
        populate_time_keywords(input_model, output_model)
        output_model.int_times = input_model.int_times.copy()

    if soss_kwargs['wave_grid_out'] is not None:
        wave_grid_model = SossWaveGridModel(wavegrid=wave_grid)
        log.info(f"Saving soss_wave_grid to {soss_kwargs['wave_grid_out']}")
        wave_grid_model.save(path=soss_kwargs['wave_grid_out'])
        wave_grid_model.close()

    return output_model, output_references, output_atoca
