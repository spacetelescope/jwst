import logging

import numpy as np
from scipy.interpolate import UnivariateSpline

from ... import datamodels
from ...datamodels import dqflags
from astropy.nddata.bitmask import bitfield_to_boolean_mask

from .soss_syscor import make_background_mask, soss_background
from .soss_solver import solve_transform, transform_wavemap, transform_profile, transform_coords
from .atoca import ExtractionEngine
from .atoca_utils import ThroughputSOSS, WebbKernel, grid_from_map, mask_bad_dispersion_direction
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


def estim_flux_first_order(scidata_bkg, scierr, scimask, ref_files, threshold):
    """
    Parameters
    ----------
    scidata_bkg : array
        A single background subtracted NIRISS SOSS detector image.
    scierr : array
        The uncertainties corresponding to the detector image.
    scimask : array
        Pixel mask to apply to the detector image.
    ref_files : list
        A list of list of the reference files for each order.
    threshold : float
        The threshold value for using pixels based on the spectral profile.

    Returns
    -------
    func
        A spline estimator that provides the underlying flux as a function of wavelength
    """

    # Unpack ref_files
    wave_maps, spat_pros, thrpts, _ = ref_files

    # Oversampling of 1 to make sure the solution will be stable
    n_os = 1

    # Define wavelength grid based on order 1 only (so first index)
    wave_grid = grid_from_map(wave_maps[0], spat_pros[0], n_os=n_os)

    # Mask parts contaminated by order 2 based on its spatial profile
    mask = (spat_pros[1] >= threshold) | scimask

    # Init extraction without convolution kernel (so extract the spectrum at order 1 resolution)
    ref_file_args = [wave_maps[0]], [spat_pros[0]], [thrpts[0]], [np.array([1.])]
    kwargs = {'wave_grid': wave_grid,
              'orders': [1],
              'global_mask': mask,
              'threshold': threshold}
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


def model_image(scidata_bkg, scierr, scimask, refmask, ref_file_args, transform=None,
                tikfac=None, n_os=5, threshold=1e-4, soss_filter='CLEAR'):
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
    ref_file_args : tuple
        A tuple of reference file arguments constructed by get_ref_file_args().
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
    soss_filter : str, optional
        Filter setting for NIRISS SOSS during observation. For use once F277W
        algorithm is determined, currently defaults to 'CLEAR'.

    Returns
    -------
    tracemodels : dict
        Dictionary of the modeled detector images for each order.
    tikfac : float
        Optimal Tikhonov factor used in extraction
    logl : float
        Log likelihood value associated with the Tikhonov factor selected.
    """

    # Some error values are 0, we need to mask those pixels for the extraction engine.
    scimask = scimask | ~(scierr > 0)

    log.info('Extracting using transformation parameters {}'.format(transform))

    # Set the c_kwargs using the minimum value of the kernels
    c_kwargs = [{'thresh': webb_ker.min_value} for webb_ker in ref_file_args[3]]

    # Initialize the Engine.
    engine = ExtractionEngine(*ref_file_args, n_os=n_os, threshold=threshold, c_kwargs=c_kwargs)

    if tikfac is None:

        log.info('Solving for the optimal Tikhonov factor.')

        # Need a rough estimate of the underlying flux to estimate the tikhonov factor
        try:
            estimate = estim_flux_first_order(scidata_bkg, scierr, scimask, ref_file_args, threshold)
            # Initial pass 8 orders of magnitude with 10 grid points.
            factors = engine.estimate_tikho_factors(estimate, log_range=[-4, 4], n_points=10)
        except Exception as e:
            log.warning(f"Error caught in first order flux estimation: {e}\n"
                        f"Using pre-defined array for flux factor estimate.")
            factors = np.logspace(-19, -10, 10)
        # Find the tikhonov factor.
        tiktests = engine.get_tikho_tests(factors, data=scidata_bkg, error=scierr, mask=scimask)
        tikfac, mode, _ = engine.best_tikho_factor(tests=tiktests, fit_mode='chi2')

        # Refine across 4 orders of magnitude.
        tikfac = np.log10(tikfac)
        factors = np.logspace(tikfac - 2, tikfac + 2, 20)
        tiktests = engine.get_tikho_tests(factors, data=scidata_bkg, error=scierr, mask=scimask)
        tikfac, mode, _ = engine.best_tikho_factor(tests=tiktests, fit_mode=mode)

    log.info('Using a Tikhonov factor of {}'.format(tikfac))

    # Run the extract method of the Engine.
    f_k = engine.__call__(data=scidata_bkg, error=scierr, mask=scimask, tikhonov=True, factor=tikfac)

    # Compute the log-likelihood of the best fit.
    logl = engine.compute_likelihood(f_k, same=False)

    log.info('Optimal solution has a log-likelihood of {}'.format(logl))

    # Create a new instance of the engine for evaluating the trace model.
    # This allows bad pixels and pixels below the threshold to be reconstructed as well.
    # Model the order 1 and order 2 trace separately.
    order_list = ['Order 1', 'Order 2']
    tracemodels = dict()
    for i_order, order in enumerate(order_list):

        log.debug('Building the model image of {}.'.format(order))

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

        # Build model of the order
        model_kwargs = {'wave_grid': grid_order,
                        'threshold': 1e-5,
                        'global_mask': refmask,
                        'orders': [i_order + 1]}
        model = ExtractionEngine(*ref_file_order, **model_kwargs)

        # Project on detector and save in dictionary
        tracemodels[order] = model.rebuild(flux_order)

    return tracemodels, tikfac, logl


def extract_image(scidata_bkg, scierr, scimask, tracemodels, ref_files,
                  transform, subarray, width=40., bad_pix='model'):
    """Perform the box-extraction on the image, while using the trace model to
    correct for contamination.
    Parameters
    ----------
    scidata_bkg : array[float]
        A single backround subtracted NIRISS SOSS detector image.
    scierr : array[float]
        The uncertainties corresponding to the detector image.
    scimask : array[float]
        Pixel mask to apply to the detector image.
    tracemodels : dict
        Dictionary of the modeled detector images for each order.
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

    Returns
    -------
    wavelengths, fluxes, fluxerrs, npixels, box_weights : dict
        Each output is a dictionary, with each extracted order as a key.
    """
    # Which orders to extract.
    if subarray == 'SUBSTRIP96':
        order_list = [1]
    else:
        order_list = [1, 2, 3]

    order_str = {order: f'Order {order}' for order in order_list}

    # List of modeled orders
    mod_order_list = tracemodels.keys()

    # Create dictionaries for the output spectra.
    wavelengths = dict()
    fluxes = dict()
    fluxerrs = dict()
    npixels = dict()
    box_weights = dict()

    log.info('Performing the decontaminated box extraction.')

    # Extract each order from order list
    for order_integer in order_list:
        # Order string-name is used more often than integer-name
        order = order_str[order_integer]

        log.debug(f'Extracting {order}.')

        # Define the box aperture
        xtrace, ytrace, wavelengths[order] = get_trace_1d(ref_files, transform, order_integer)
        box_w_ord = get_box_weights(ytrace, width, scidata_bkg.shape, cols=xtrace)

        # Decontaminate using all other modeled orders
        decont = scidata_bkg
        for mod_order in mod_order_list:
            if mod_order != order:
                log.debug(f'Decontaminating {order} from {mod_order} using model.')
                decont = decont - tracemodels[mod_order]

        # Deal with bad pixels if required.
        if bad_pix == 'model':
            # Model the bad pixels decontaminated image when available
            try:
                # Replace bad pixels
                decont = np.where(scimask, tracemodels[order], decont)
                # Update the mask for the modeled order, so all the pixels are usable.
                scimask_ord = np.zeros_like(scimask)

                log.debug(f'Bad pixels in {order} are replaced with trace model.')

                # Replace error estimate of the bad pixels using other valid pixels of similar value.
                # The pixel to be estimate are the masked pixels in the region of extraction
                extraction_region = (box_w_ord > 0)
                pix_to_estim = (extraction_region & scimask)
                # Use only valid pixels (not masked) in the extraction region for the empirical estimation
                valid_pix = (extraction_region & ~scimask)
                scierr_ord = estim_error_nearest_data(scierr, decont, pix_to_estim, valid_pix)

            except KeyError:
                # Keep same mask and error
                scimask_ord = scimask
                scierr_ord = scierr
                log.warning(f'Bad pixels in {order} will be masked instead of modeled: trace model unavailable.')
        else:
            # Mask pixels
            scimask_ord = scimask
            scierr_ord = scierr
            log.info(f'Bad pixels in {order} will be masked.')

        # Save box weights
        box_weights[order] = box_w_ord
        # Perform the box extraction and save
        out = box_extract(decont, scierr_ord, scimask_ord, box_w_ord, cols=xtrace)
        _, fluxes[order], fluxerrs[order], npixels[order] = out

    return wavelengths, fluxes, fluxerrs, npixels, box_weights


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

    # Initialize the output model and output references (model of the detector and box aperture weights).
    output_model = datamodels.MultiSpecModel()
    output_model.update(input_model)  # Copy meta data from input to output.

    output_references = datamodels.SossExtractModel()
    output_references.update(input_model)

    all_tracemodels = dict()
    all_box_weights = dict()

    # Extract depending on the type of datamodels (Image or Cube)
    if isinstance(input_model, datamodels.ImageModel):

        log.info('Input is an ImageModel, processing a single integration.')

        # Initialize the theta, dx, dy transform parameters
        transform = soss_kwargs.pop('transform')

        # Received a single 2D image; set dtype to float64 and convert DQ to boolean mask.
        scidata = input_model.data.astype('float64')
        scierr = input_model.err.astype('float64')
        scimask = input_model.dq > 0  # Mask bad pixels with True.
        refmask = bitfield_to_boolean_mask(input_model.dq,
                                           ignore_flags=dqflags.pixel['REFERENCE_PIXEL'],
                                           flip_bits=True)

        # Perform background correction.
        bkg_mask = make_background_mask(scidata, width=40)
        scidata_bkg, col_bkg, npix_bkg = soss_background(scidata, scimask, bkg_mask=bkg_mask)

        # Determine the theta, dx, dy transform needed to match scidata trace position to ref file position.
        if transform is None:
            log.info('Solving for the transformation parameters.')

            # Unpack the expected order 1 & 2 positions.
            spectrace_ref = ref_files['spectrace']
            xref_o1 = spectrace_ref.trace[0].data['X']
            yref_o1 = spectrace_ref.trace[0].data['Y']
            xref_o2 = spectrace_ref.trace[1].data['X']
            yref_o2 = spectrace_ref.trace[1].data['Y']

            # Use the solver on the background subtracted image.
            if subarray == 'SUBSTRIP96' or soss_filter == 'F277W':
                # Use only order 1 to solve theta, dx, dy
                transform = solve_transform(scidata_bkg, scimask, xref_o1, yref_o1, soss_filter=soss_filter)
            else:
                transform = solve_transform(scidata_bkg, scimask, xref_o1, yref_o1,
                                            xref_o2, yref_o2, soss_filter=soss_filter)

        log.info('Measured to Reference trace position transform: theta={:.4f}, dx={:.4f}, dy={:.4f}'.format(
            transform[0], transform[1], transform[2]))

        # Prepare the reference file arguments.
        ref_file_args = get_ref_file_args(ref_files, transform)

        # Make sure wavelength maps cover only parts where the centroid is inside the detector image
        if subarray != 'SUBSTRIP96':
            _mask_wv_map_centroid_outside(ref_file_args[0], ref_files, transform, scidata_bkg.shape[0])

        # Model the traces based on optics filter configuration (CLEAR or F277W)
        if soss_filter == 'CLEAR':

            # Model the image.
            kwargs = dict()
            kwargs['transform'] = transform
            kwargs['tikfac'] = soss_kwargs['tikfac']
            kwargs['n_os'] = soss_kwargs['n_os']
            kwargs['threshold'] = soss_kwargs['threshold']

            result = model_image(scidata_bkg, scierr, scimask, refmask, ref_file_args, **kwargs)
            tracemodels, soss_kwargs['tikfac'], logl = result

        else:
            # No model can be fit for F277W yet, missing throughput reference files.
            msg = f"No extraction possible for filter {soss_filter}."
            log.critical(msg)
            raise ValueError(msg)

        # Save trace models for output reference
        for order in tracemodels:
            # Save as a list (convert to array at the end)
            all_tracemodels[order] = [tracemodels[order]]

        # Use the trace models to perform a decontaminated extraction.
        kwargs = dict()
        kwargs['width'] = soss_kwargs['width']
        kwargs['bad_pix'] = soss_kwargs['bad_pix']

        result = extract_image(scidata_bkg, scierr, scimask, tracemodels, ref_files, transform, subarray, **kwargs)
        wavelengths, fluxes, fluxerrs, npixels, box_weights = result

        # Save box weights for output reference
        for order in box_weights:
            # Save as a list (convert to array at the end)
            all_box_weights[order] = [box_weights[order]]

        # Copy spectral data for each order into the output model.
        for order in wavelengths.keys():

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

            output_model.spec.append(spec)

        output_model.meta.soss_extract1d.width = kwargs['width']
        output_model.meta.soss_extract1d.tikhonov_factor = soss_kwargs['tikfac']
        output_model.meta.soss_extract1d.delta_x = transform[1]
        output_model.meta.soss_extract1d.delta_y = transform[2]
        output_model.meta.soss_extract1d.theta = transform[0]
        output_model.meta.soss_extract1d.oversampling = soss_kwargs['n_os']
        output_model.meta.soss_extract1d.threshold = soss_kwargs['threshold']

    elif isinstance(input_model, datamodels.CubeModel):

        nimages = len(input_model.data)

        log.info('Input is a CubeModel containing {} integrations.'.format(nimages))

        # Initialize the theta, dx, dy transform parameters
        transform = soss_kwargs.pop('transform')

        # Loop over images.
        for i in range(nimages):

            log.info('Processing integration {} of {}.'.format(i + 1, nimages))

            # Unpack the i-th image, set dtype to float64 and convert DQ to boolean mask.
            scidata = input_model.data[i].astype('float64')
            scierr = input_model.err[i].astype('float64')
            scimask = np.bitwise_and(input_model.dq[i], dqflags.pixel['DO_NOT_USE']).astype(bool)
            refmask = bitfield_to_boolean_mask(input_model.dq[i], ignore_flags=dqflags.pixel['REFERENCE_PIXEL'],
                                               flip_bits=True)

            # Perform background correction.
            bkg_mask = make_background_mask(scidata, width=40)
            scidata_bkg, col_bkg, npix_bkg = soss_background(scidata, scimask, bkg_mask=bkg_mask)

            # Determine the theta, dx, dy transform needed to match scidata trace position to ref file position.
            if transform is None:
                log.info('Solving for the transformation parameters.')

                # Unpack the expected order 1 & 2 positions.
                spectrace_ref = ref_files['spectrace']
                xref_o1 = spectrace_ref.trace[0].data['X']
                yref_o1 = spectrace_ref.trace[0].data['Y']
                xref_o2 = spectrace_ref.trace[1].data['X']
                yref_o2 = spectrace_ref.trace[1].data['Y']

                # Use the solver on the background subtracted image.
                if subarray == 'SUBSTRIP96' or soss_filter == 'F277W':
                    # Use only order 1 to solve theta, dx, dy
                    transform = solve_transform(scidata_bkg, scimask, xref_o1, yref_o1, soss_filter=soss_filter)
                else:
                    transform = solve_transform(scidata_bkg, scimask, xref_o1, yref_o1,
                                                xref_o2, yref_o2, soss_filter=soss_filter)

            log.info('Measured to Reference trace position transform: theta={:.4f}, dx={:.4f}, dy={:.4f}'.format(
                     transform[0], transform[1], transform[2]))

            # Prepare the reference file arguments.
            ref_file_args = get_ref_file_args(ref_files, transform)

            # Make sure wavelength maps cover only parts where the centroid is inside the detector image
            _mask_wv_map_centroid_outside(ref_file_args[0], ref_files, transform, scidata_bkg.shape[0])

            # Model the traces based on optics filter configuration (CLEAR or F277W)
            if soss_filter == 'CLEAR':

                # Model the image.
                kwargs = dict()
                kwargs['transform'] = transform
                kwargs['tikfac'] = soss_kwargs['tikfac']
                kwargs['n_os'] = soss_kwargs['n_os']
                kwargs['threshold'] = soss_kwargs['threshold']

                result = model_image(scidata_bkg, scierr, scimask, refmask, ref_file_args, **kwargs)
                tracemodels, soss_kwargs['tikfac'], logl = result

            else:
                # No model can be fit for F277W yet, missing throughput reference files.
                msg = f"No extraction possible for filter {soss_filter}."
                log.critical(msg)
                return None, None

            # Save trace models for output reference
            for order in tracemodels:
                # Initialize a list for first integration
                if i == 0:
                    all_tracemodels[order] = []
                all_tracemodels[order].append(tracemodels[order])

            # Use the trace models to perform a de-contaminated extraction.
            kwargs = dict()
            kwargs['width'] = soss_kwargs['width']
            kwargs['bad_pix'] = soss_kwargs['bad_pix']

            result = extract_image(scidata_bkg, scierr, scimask, tracemodels, ref_files, transform, subarray, **kwargs)
            wavelengths, fluxes, fluxerrs, npixels, box_weights = result

            # Save box weights for output reference
            for order in box_weights:
                # Initialize a list for first integration
                if i == 0:
                    all_box_weights[order] = []
                all_box_weights[order].append(box_weights[order])

            # Copy spectral data for each order into the output model.
            for order in wavelengths.keys():

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

            output_model.meta.soss_extract1d.width = kwargs['width']
            output_model.meta.soss_extract1d.tikhonov_factor = soss_kwargs['tikfac']
            output_model.meta.soss_extract1d.delta_x = transform[1]
            output_model.meta.soss_extract1d.delta_y = transform[2]
            output_model.meta.soss_extract1d.theta = transform[0]
            output_model.meta.soss_extract1d.oversampling = soss_kwargs['n_os']
            output_model.meta.soss_extract1d.threshold = soss_kwargs['threshold']

    else:
        msg = "Only ImageModel and CubeModel are implemented for the NIRISS SOSS extraction."
        log.critical(msg)
        raise ValueError(msg)

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

    return output_model, output_references
