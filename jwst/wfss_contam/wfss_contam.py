"""Top-level module for WFSS contamination correction."""

import logging
import multiprocessing
import types

import numpy as np
from stcal.multiprocessing import compute_num_cores
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.transforms.models import (
    NIRCAMBackwardGrismDispersion,
    NIRISSBackwardGrismDispersion,
)

from jwst.lib.catalog_utils import read_source_catalog
from jwst.wfss_contam.observations import Observation
from jwst.wfss_contam.sens1d import get_photom_data
from jwst.wfss_contam.wavefit import SlitFitError, apply_basis_coeffs, fit_slit_by_basis_images

log = logging.getLogger(__name__)

__all__ = ["contam_corr"]


class UnmatchedSlitIDError(Exception):
    """Exception raised when a slit ID is not found in the list of simulated slits."""

    pass


class _LegendreFluxModel:
    """
    Picklable callable that evaluates the k-th Legendre polynomial.

    Callables need to be picklable so that they can be used in multiprocessing contexts.
    The wavelength argument is mapped from ``[wmin, wmax]`` to ``[-1, 1]`` before
    evaluation.
    """

    def __init__(self, k, wmin, wmax):
        self.k = k
        self.wmin = wmin
        self.wmax = wmax
        # Coefficient vector with 1 at position k and 0 elsewhere.
        self._coeffs = np.zeros(k + 1)
        self._coeffs[k] = 1.0

    def __call__(self, x):
        x_norm = 2.0 * (x - self.wmin) / (self.wmax - self.wmin) - 1.0
        return np.polynomial.legendre.legval(x_norm, self._coeffs)


def _find_matching_simul_slit(slit, simul_slit_sids, simul_slit_orders):
    """
    Find the index of the matching simulated slit in the list of simulated slits.

    Parameters
    ----------
    slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Source slit model
    simul_slit_sids : list
        List of source IDs for simulated slits
    simul_slit_orders : list
        List of spectral orders for simulated slits

    Returns
    -------
    good_idx : int
        Index of the matching simulated slit in the list of simulated slits
    """
    sid = slit.source_id
    order = slit.meta.wcsinfo.spectral_order
    good = (np.array(simul_slit_sids) == sid) * (np.array(simul_slit_orders) == order)
    if not any(good):
        raise UnmatchedSlitIDError(
            f"Source ID {sid} order {order} requested by input slit model "
            "but not found in simulated slits. "
            "Setting contamination correction to zero for that slit."
        )
    return np.where(good)[0][0]


def _cut_frame_to_match_slit(contam, slit):
    """
    Cut out the contamination image to match the extent of the source slit.

    Parameters
    ----------
    contam : 2D array
        Contamination image for the full grism exposure
    slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Source slit model

    Returns
    -------
    cutout : 2D array
        Contamination image cutout that matches the extent of the source slit
    """
    x1 = slit.xstart - 1  # convert FITS 1-indexed to 0-indexed for array access
    y1 = slit.ystart - 1  # convert FITS 1-indexed to 0-indexed for array access
    xf = x1 + slit.xsize
    yf = y1 + slit.ysize

    # zero-pad the contamination image if the slit extends beyond the contamination image
    # fixes an off-by-one bug when sources extend to the edge of the contamination image
    if xf > contam.shape[1]:
        contam = np.pad(contam, ((0, 0), (0, xf - contam.shape[1])), mode="constant")
    if yf > contam.shape[0]:
        contam = np.pad(contam, ((0, yf - contam.shape[0]), (0, 0)), mode="constant")

    return contam[y1 : y1 + slit.ysize, x1 : x1 + slit.xsize]


class SlitOverlapError(Exception):
    """Exception raised when there is no overlap between data and model for a slit."""

    pass


def match_backplane_prefer_first(slit0, slit1):
    """
    Reshape slit1 to the backplane of slit0.

    Parameters
    ----------
    slit0 : `~stdatamodels.jwst.datamodels.SlitModel`
        Slit model for the first slit, which is used as reference.
    slit1 : `~stdatamodels.jwst.datamodels.SlitModel`
        Slit model for the second slit, which is reshaped to match slit0.

    Returns
    -------
    slit1 : `~stdatamodels.jwst.datamodels.SlitModel`
        Reshaped slit model.
    """
    data0 = slit0.data
    data1 = slit1.data

    x1 = slit1.xstart - slit0.xstart
    y1 = slit1.ystart - slit0.ystart
    backplane1 = np.zeros_like(data0)

    i0 = max([y1, 0])
    i1 = min([y1 + data1.shape[0], data0.shape[0], data1.shape[0]])
    j0 = max([x1, 0])
    j1 = min([x1 + data1.shape[1], data0.shape[1], data1.shape[1]])
    if i0 >= i1 or j0 >= j1:
        raise SlitOverlapError(
            f"No overlap region between data and model for slit {slit0.source_id}, "
            f"order {slit0.meta.wcsinfo.spectral_order}. "
            "setting contamination correction to zero for that slit."
        )
    di = i0 - y1  # offset into data1's own row axis
    dj = j0 - x1  # offset into data1's own col axis
    backplane1[i0:i1, j0:j1] = data1[di : di + (i1 - i0), dj : dj + (j1 - j0)]
    slit1.data = backplane1

    # also update fluxmodel attributes if present
    k = 1
    while True:
        attr = f"fluxmodel_{k}"
        mc = getattr(slit1, attr, None)
        if mc is None:
            break
        mc_bp = np.zeros_like(data0)
        mc_bp[i0:i1, j0:j1] = np.asarray(mc)[di : di + (i1 - i0), dj : dj + (j1 - j0)]
        setattr(slit1, attr, mc_bp)
        k += 1

    slit1.xstart = slit0.xstart
    slit1.ystart = slit0.ystart
    slit1.xsize = slit0.xsize
    slit1.ysize = slit0.ysize

    return slit1


def _validate_orders_against_reference(orders, spec_orders):
    """
    Compare user-requested spectral orders with the orders defined in the reference file.

    Parameters
    ----------
    orders : list[int]
        List of user-requested spectral orders.
    spec_orders : list[int]
        List of spectral orders defined in the reference file.

    Returns
    -------
    np.ndarray[int]
        List of spectral orders constrained to the user-specified ones
        that are also defined in the reference file.
    """
    spec_orders = np.array(spec_orders, dtype=int)
    if orders is None:
        return spec_orders
    orders = np.array(orders, dtype=int)
    good_orders = np.isin(orders, spec_orders, assume_unique=True)
    if (len(good_orders) == 0) or (not np.any(good_orders)):
        log.error(
            f"None of the requested spectral orders {orders} are defined "
            "in the wavelength range reference file. "
            f"Expected orders are: {spec_orders}. "
        )
        return []
    if not np.all(good_orders):
        log.warning(
            f"Not all requested spectral orders {orders} are defined in the "
            f"wavelength range reference file. Defined orders are: {spec_orders}. "
            "Skipping undefined orders."
        )
    return orders[good_orders]


def _validate_orders_against_transform(wcs, spec_orders):
    """
    Ensure the requested spectral orders are defined in the WCS transforms.

    Parameters
    ----------
    wcs : gwcs.wcs.WCS
        The input MultiSlitModel's WCS object.
    spec_orders : list[int]
        The list of requested spectral orders.

    Returns
    -------
    list
        List of spectral orders that are defined in the WCS transform.
    """
    sky_to_grism = wcs.backward_transform
    good_orders = spec_orders.copy()
    for model in sky_to_grism:
        if isinstance(model, (NIRCAMBackwardGrismDispersion, NIRISSBackwardGrismDispersion)):
            # Get the orders defined in the transform
            orders = np.sort(model.orders)
            is_good_order = [order in orders for order in spec_orders]
            if not any(is_good_order):
                log.error(
                    f"None of the requested spectral orders {spec_orders} are defined "
                    "in the WCS transform. "
                    f"Defined orders are: {orders}. "
                )
                return []
            if not all(is_good_order):
                log.warning(
                    f"Not all requested spectral orders {spec_orders} are "
                    f"defined in the WCS transform. Defined orders are: {orders}. "
                    "Skipping undefined orders."
                )
            good_orders = [order for order in spec_orders if order in orders]
            # There will be only one transform of this type in the wcs
            break
    return np.sort(good_orders)


def _find_min_relresp(sens_waves, sens_response):
    """
    Find the minimum relative response in the sensitivity response.

    Helper function is necessary instead of just nanmin because reference file
    sometimes has zero-valued wavelength/response pairs.

    Parameters
    ----------
    sens_waves : np.ndarray
        Wavelengths corresponding to the sensitivity response.
    sens_response : np.ndarray
        Sensitivity response values.

    Returns
    -------
    float
        Minimum relative response value.
    """
    good = (sens_waves > 0) & (sens_response > 0) & np.isfinite(sens_waves)
    return np.nanmin(sens_response[good])


def _build_simulated_image_from_slits(simulated_slits, shape):
    """
    Reconstruct the full-frame simulated image from the simulated slits.

    This is needed instead of just using ``obs.simulated_image`` because when
    spectral fitting is requested, the simulated slits get modified, and the
    modifications need to end up in the simul file.

    Parameters
    ----------
    simulated_slits : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        The simulated slits.
    shape : tuple of int
        ``(nrows, ncols)`` of the full detector frame.

    Returns
    -------
    full_image : ndarray
        Full-frame simulated image.
    """
    full_image = np.zeros(shape, dtype=float)
    nrows, ncols = shape
    for slit in simulated_slits.slits:
        x0 = slit.xstart - 1  # convert FITS 1-indexed to 0-indexed for array access
        y0 = slit.ystart - 1  # convert FITS 1-indexed to 0-indexed for array access
        # Clip to frame boundaries in case a slit overflows
        x1 = min(x0 + slit.xsize, ncols)
        y1 = min(y0 + slit.ysize, nrows)
        full_image[y0:y1, x0:x1] += slit.data[: y1 - y0, : x1 - x0]
    return full_image


def _apply_magnitude_limit(
    order, source_catalog, sens_wave, sens_response, magnitude_limit, min_relresp_order1
):
    """
    Rescale the magnitude limit based on the sensitivity response for a given spectral order.

    Parameters
    ----------
    order : int
        Spectral order for which the magnitude limit is applied.
    source_catalog : astropy.table.Table
        The source catalog containing source IDs and isophotal AB magnitudes.
    sens_wave : np.ndarray
        The wavelengths corresponding to the sensitivity response.
    sens_response : np.ndarray
        The sensitivity response for the order.
    magnitude_limit : float
        The isophotal AB magnitude limit for sources to be included in the contamination correction.
    min_relresp_order1 : float
        Minimum relative response for order 1, used to scale the magnitude limit.

    Returns
    -------
    list
        List of source IDs that meet the magnitude limit criteria.
    """
    if order in [0, 1]:
        # Magnitude limit is set according to order 1 sensitivity response
        # and order 0 is a special case because it's not dispersed
        order_mag_limit = magnitude_limit
    else:
        # Scale the magnitude limit according to the order sensitivity response
        order_sens_factor = min_relresp_order1 / _find_min_relresp(sens_wave, sens_response)
        order_mag_diff = -2.5 * np.log10(order_sens_factor)
        order_mag_limit = magnitude_limit - order_mag_diff

    # Select sources that are brighter than the magnitude limit
    good_sources = source_catalog[source_catalog["isophotal_abmag"] < order_mag_limit]
    if len(good_sources) == 0:
        return None
    log.info(
        f"Applying magnitude limit of {order_mag_limit:.1f} to order {order}. "
        f"Sources selected: {len(good_sources)}"
    )
    return good_sources["label"].tolist()


def _match_simulated_slits(output_model, obs):
    """
    Match each observed slit to its simulated counterpart.

    Reprojects the simulated slit onto the same backplane as the observed slit.
    If the slit ID is not found in the simulated slits, or if there is no spatial overlap between
    the observed slit and the simulated slit, the slit is skipped and represented as ``None``
    in the returned list.

    Parameters
    ----------
    output_model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        The observed slits to match against.
    obs : `~jwst.wfss_contam.observations.Observation`
        The observation object containing the dispersed simulated slits.

    Returns
    -------
    matched_flat_simuls : list of `~stdatamodels.jwst.datamodels.SlitModel`
        Simulated slit reprojected onto each observed slit's backplane,
        or ``None`` where no match was found.
    good_idxs : list of int or None
        Index into ``obs.simulated_slits.slits`` for each observed slit,
        or ``None`` where no match was found.
    """
    simul_slit_sids = [slit.source_id for slit in obs.simulated_slits.slits]
    simul_slit_orders = [slit.meta.wcsinfo.spectral_order for slit in obs.simulated_slits.slits]

    matched_flat_simuls = []
    good_idxs = []
    for slit in output_model.slits:
        try:
            good_idx = _find_matching_simul_slit(slit, simul_slit_sids, simul_slit_orders)
            # Copy only the data arrays and metadata specifically needed by
            # match_backplane_prefer_first and fit_slit_by_basis_images.
            # This reduces overall peak memory usage of the step by a factor of ~3 in some cases
            src = obs.simulated_slits.slits[good_idx]
            matched_flat = types.SimpleNamespace(
                data=np.array(src.data),
                xstart=src.xstart,
                ystart=src.ystart,
                xsize=src.xsize,
                ysize=src.ysize,
            )
            k = 1
            while True:
                mc = getattr(src, f"fluxmodel_{k}", None)
                if mc is None:
                    break
                setattr(matched_flat, f"fluxmodel_{k}", np.array(mc))
                k += 1
            matched_flat = match_backplane_prefer_first(slit, matched_flat)
            matched_flat_simuls.append(matched_flat)
            good_idxs.append(good_idx)
        except (UnmatchedSlitIDError, SlitOverlapError) as e:
            log.warning(e)
            matched_flat_simuls.append(None)
            good_idxs.append(None)

    return matched_flat_simuls, good_idxs


def _fit_spectral_shape(
    observed_slit,
    simul_slit,
    simul_slit_backplane_unmatched,
    polyfit_degree,
    l2_alpha=0.0,
    rejection_threshold=0.1,
):
    """
    Fit a polynomial spectral shape to one slit and apply the result in-place.

    Parameters
    ----------
    observed_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Observed slit whose ``.data`` is used as the target for the fit.
        If ``n_iterations > 1``, this holds the contamination-corrected data
        from the previous iteration.
    simul_slit : `~stdatamodels.jwst.datamodels.SlitModel`
        Simulated slit, already backplane-matched to ``observed_slit``.
        Its ``.data`` attribute is updated in-place with the spectrally fitted result.
    simul_slit_backplane_unmatched : `~stdatamodels.jwst.datamodels.SlitModel`
        The simulation in ``obs.simulated_slits`` without backplane matching applied.
        This is tracked independently from ``simul_slit`` because the extraction
        of the observed slit from extract_2d often misses flux from the extended PSF wings,
        whereas the simulation in this step disperses the whole segment, including those wings.
        Its ``.data`` is updated in-place so the full-frame reconstruction
        reflects the fitted spectral shape.
    polyfit_degree : int or None
        Degree of the polynomial spectral model.  ``None`` means no fitting.
    l2_alpha : float, optional
        L2 regularisation strength passed to `~jwst.wfss_contam.wavefit.fit_slit_by_basis_images`.
    rejection_threshold : float, optional
        Threshold for rejecting fits based on the fitted constant term coefficient, passed to
        `~jwst.wfss_contam.wavefit.fit_slit_by_basis_images`.

    Returns
    -------
    status : bool
        ``True`` if the fit was successful, ``False`` if the fit failed and the
        original flat-spectrum simulation is used.
    """
    if polyfit_degree is None or getattr(simul_slit, "fluxmodel_1", None) is None:
        return False

    log.debug(
        f"Fitting polynomial of degree {polyfit_degree} to the simulated slit "
        f"for source ID {observed_slit.source_id}, "
        f"order {observed_slit.meta.wcsinfo.spectral_order}"
    )
    try:
        coeffs = fit_slit_by_basis_images(
            observed_slit, simul_slit, l2_alpha=l2_alpha, rejection_threshold=rejection_threshold
        )
        if coeffs is None:
            return False
        simul_slit.data = apply_basis_coeffs(simul_slit, coeffs)
        simul_slit_backplane_unmatched.data = apply_basis_coeffs(
            simul_slit_backplane_unmatched, coeffs
        )
    except SlitFitError as e:
        log.debug(
            f"Polynomial fitting failed for slit with source ID {observed_slit.source_id}, "
            f"order {observed_slit.meta.wcsinfo.spectral_order}: {e}. "
            "Using the original simulated slit without fitting."
        )
        return False
    else:
        return True


def _build_contam(output_model, per_slit_simuls, simul_data, original_data):
    """
    Build the contamination model for each slit.

    Parameters
    ----------
    output_model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        The output model containing the observed spectral cutouts.
    per_slit_simuls : list
        List of simulated spectra corresponding to each observed cutout.
    simul_data : `~stdatamodels.jwst.datamodels.SlitModel`
        The full-frame simulated data.
    original_data : list
        List of original observed data arrays.

    Returns
    -------
    list
        List of contamination cutouts for each slit.
    """
    contam_cuts = []
    for i, (slit, this_simul) in enumerate(zip(output_model.slits, per_slit_simuls, strict=True)):
        if this_simul is None:
            contam_cut = np.zeros_like(slit.data)
        else:
            simul_all_cut = _cut_frame_to_match_slit(simul_data, slit)
            contam_cut = simul_all_cut - this_simul.data
        slit.data = original_data[i] - contam_cut
        contam_cuts.append(contam_cut)
    return contam_cuts


def contam_corr(
    input_model,
    waverange,
    photom,
    max_cores,
    orders=None,
    magnitude_limit=None,
    max_pixels_per_chunk=5e4,
    oversample_factor=2,
    polyfit_degree=None,
    n_iterations=1,
    l2_alpha=0.1,
    rejection_threshold=0.1,
):
    """
    Correct contamination in WFSS spectral cutouts.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        Input data model containing 2D spectral cutouts. May be modified by processing:
        make a copy before calling this function, if needed.
    waverange : `~stdatamodels.jwst.datamodels.WavelengthrangeModel`
        Wavelength range reference file model
    photom : `~stdatamodels.jwst.datamodels.NrcWfssPhotomModel` or \
             `~stdatamodels.jwst.datamodels.NisWfssPhotomModel`
        Photom (flux cal) reference file model
    max_cores : str or int
        Number of cores to use for multiprocessing. If set to 'none'
        (the default), then no multiprocessing will be done. The other
        allowable string values are 'quarter', 'half', and 'all', which indicate
        the fraction of cores to use for multi-proc. The total number of
        cores includes the SMT cores (Hyper Threading for Intel).
        If an integer is provided, it will be the exact number of cores used.
    orders : list, optional
        List of spectral orders to process.
        If None, all orders defined in the wavelengthrange file will be processed.
    magnitude_limit : float, optional
        Isophotal AB magnitude limit for sources to be included in the contamination correction.
        The magnitude limit is applied per spectral order, where the orders are scaled relative
        to order 0 based on their photometric response as read from the photom reference file.
        This means that generally fewer sources will be dispersed in higher orders.
        If None, no magnitude limit is applied and all sources are included.
    max_pixels_per_chunk : int, optional
        Maximum number of pixels to disperse simultaneously.
    oversample_factor : int, optional
        Wavelength oversampling factor.
    polyfit_degree : int, optional
        Degree of polynomial fit to spectral shape. If None (the default), do not attempt
        polynomial fitting and just use the flat-spectrum simulated slit.
    n_iterations : int, optional
        Number of times to iterate the contamination correction. On each iteration the
        polynomial fit is re-run using the contamination-corrected spectrum from the
        previous iteration, yielding a progressively better estimate of each source's
        true spectral shape (and therefore a better contamination estimate for its
        neighbors). Requires ``polyfit_degree`` to be set; if ``polyfit_degree`` is
        None this parameter is ignored and a single iteration is performed.
    l2_alpha : float, optional
        L2 regularization strength for the polynomial spectral fit, passed to
        `~jwst.wfss_contam.wavefit.fit_slit_by_basis_images`.
    rejection_threshold : float, optional
        Threshold for rejecting fits based on the fitted constant term coefficient, passed to
        `~jwst.wfss_contam.wavefit.fit_slit_by_basis_images`.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        A copy of the input_model that has been corrected
    simul_model : `~stdatamodels.jwst.datamodels.ImageModel`
        Full-frame simulated image of the grism exposure
    contam_model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
        Contamination estimate images for each source slit
    """
    max_available_cores = multiprocessing.cpu_count()
    # don't worry about case where nchunks < ncpus; just set nchunks large for now
    ncpus = compute_num_cores(max_cores, 1e10, max_available_cores)

    # Get the segmentation map and direct image for this grism exposure
    seg_model = datamodels.open(input_model.meta.segmentation_map)
    direct_file = input_model.meta.direct_image
    log.debug(f"Direct image ={direct_file}")
    with datamodels.open(direct_file) as direct_model:
        direct_image = direct_model.data
        direct_image_wcs = direct_model.meta.wcs

    # Get the grism WCS object from the first cutout in the input model.
    # This WCS is used to transform from direct image to grism frame for all sources
    # in the segmentation map.
    # The "detector" to "grism_detector" and "world" to "detector" transforms are identical
    # for all slits, so just use the first one. The "grism_detector" to "grism_slit"
    # transform is not used by the step.
    grism_wcs = input_model.slits[0].meta.wcs

    # Find out how many spectral orders are defined based on the
    # array of order values in the Wavelengthrange ref file,
    # then constrain the orders to the user-specified ones
    spec_orders = np.asarray(waverange.order)
    spec_orders = _validate_orders_against_reference(orders, spec_orders)
    spec_orders = _validate_orders_against_transform(grism_wcs, spec_orders)
    if len(spec_orders) == 0:
        log.error("No valid spectral orders found. Step will be SKIPPED.")
        return input_model, None, None, None
    log.info(f"Spectral orders requested = {[int(x) for x in spec_orders]}")

    # Get the FILTER and PUPIL wheel positions, for use later
    filter_kwd = input_model.meta.instrument.filter
    pupil_kwd = input_model.meta.instrument.pupil

    # NOTE: The NIRCam WFSS mode uses filters that are in the FILTER wheel
    # with gratings in the PUPIL wheel. NIRISS WFSS mode, however, is just
    # the opposite. It has gratings in the FILTER wheel and filters in the
    # PUPIL wheel. So when processing NIRISS grism exposures the name of
    # filter needs to come from the PUPIL keyword value.
    if input_model.meta.instrument.name == "NIRISS":
        filter_name = pupil_kwd
    else:
        filter_name = filter_kwd

    # Read the source catalog to perform magnitude-based source selection later
    # mag limit will be scaled according to order 1 sensitivity
    if magnitude_limit is not None:
        source_catalog = read_source_catalog(input_model.meta.source_catalog)
        order1_wave_response, order1_sens_response = get_photom_data(
            photom, filter_kwd, pupil_kwd, order=1
        )
        min_relresp_order1 = _find_min_relresp(order1_wave_response, order1_sens_response)

    # set up observation object to disperse
    obs = Observation(
        direct_image,
        seg_model.data,
        grism_wcs,
        direct_image_wcs,
        boundaries=[0, 2047, 0, 2047],
        max_cpu=ncpus,
        max_pixels_per_chunk=max_pixels_per_chunk,
        oversample_factor=oversample_factor,
    )
    seg_model.close()

    no_sources = True
    for order in spec_orders:
        # Load lists of wavelength ranges and flux cal info
        wavelength_range = waverange.get_wfss_wavelength_range(filter_name, [order])
        wmin = wavelength_range[order][0]
        wmax = wavelength_range[order][1]
        log.debug(f"wmin={wmin}, wmax={wmax} for order {order}")
        sens_waves, sens_response = get_photom_data(photom, filter_kwd, pupil_kwd, order)

        # Build Legendre basis flux models for disperse() if polynomial fitting is requested.
        # wmin/wmax are order-specific, so construction must happen inside the order loop.
        # The constant term (k=0, i.e. slit.data) is always included; start at degree 1.
        basis_models = None
        if polyfit_degree is not None:
            basis_models = [_LegendreFluxModel(k, wmin, wmax) for k in range(1, polyfit_degree + 1)]

        # Constrain the source IDs to those that are below the magnitude limit
        selected_ids = None
        if magnitude_limit is not None:
            good_ids = _apply_magnitude_limit(
                order,
                source_catalog,
                sens_waves,
                sens_response,
                magnitude_limit,
                min_relresp_order1,
            )
            if good_ids is None:
                log.info(
                    f"No sources meet the magnitude limit of {magnitude_limit} for order {order}. "
                    "Skipping contamination correction for this order."
                )
                continue
            selected_ids = good_ids
        no_sources = False

        # Compute the dispersion for all sources in this order
        log.info(f"Creating full simulated grism image for order {order}")
        obs.disperse_order(
            order,
            wmin,
            wmax,
            sens_waves,
            sens_response,
            selected_ids,
            basis_models=basis_models,
        )

    if no_sources:
        log.error(
            f"No sources found that met the magnitude limit {magnitude_limit}. Step will be SKIPPED"
        )
        return input_model, None, None, None

    # Initialize output multislitmodel
    output_model = datamodels.MultiSlitModel()

    # Copy over matching slits.
    # Note that this makes a reference to input slits, not a deep copy,
    # so the input data may be modified by this function.  The input data is
    # copied in the calling step, as needed.
    good_slits = [slit for slit in input_model.slits if slit.source_id in obs.source_ids]
    output_model.slits.extend(good_slits)

    contam_model = datamodels.MultiSlitModel()
    contam_model.update(input_model, only="PRIMARY")
    simul_slits = datamodels.MultiSlitModel()
    simul_slits.update(input_model, only="PRIMARY")

    # Hold onto original input data so iterative corrections always start from the same baseline.
    original_data = [np.array(slit.data) for slit in output_model.slits]

    if n_iterations > 1 and polyfit_degree is None:
        log.warning(
            "n_iterations > 1 has no effect when polyfit_degree is None "
            "(there is no spectral fit to iterate). Only one iteration will be performed."
        )
        n_iterations = 1

    # Match simulated slits to observed slits
    matched_flat_simuls, good_idxs = _match_simulated_slits(output_model, obs)

    if polyfit_degree is not None:
        # Iterate: each pass re-fits spectral shapes using the contamination-corrected
        # data from the previous pass, giving progressively better contamination estimates.
        log.info(
            f"Using polyfit_degree={polyfit_degree} "
            f"for spectral fitting over {n_iterations} iterations"
        )
        # check that background subtraction did not fail
        if input_model.meta.cal_step.bkg_subtract != "COMPLETE":
            log.warning(
                f"Background subtraction step status is {input_model.meta.cal_step.bkg_subtract}. "
                "A good background subtraction is necessary for models to match observed data "
                "well enough for spectral fitting to succeed. Fitting will be attempted, "
                "but failures may be expected."
            )

    # Save the brightness of the simulation for each spectrum for sorting later.
    # _fit_spectral_shape only reassigns .data (never modifies fluxmodel_N), so
    # the flat data captures the intrinsic source flux independently of any fitted shape.
    flat_matched_sum = [
        np.nansum(np.array(s.data)) if s is not None else None for s in matched_flat_simuls
    ]

    # Apply flat-spectrum contamination correction
    # If fitting is requested, fit will start with flat-contam-removed data
    per_slit_simuls = list(matched_flat_simuls)
    simul_data = _build_simulated_image_from_slits(obs.simulated_slits, obs.simulated_image.shape)
    contam_cuts = _build_contam(output_model, per_slit_simuls, simul_data, original_data)

    if polyfit_degree is not None:
        for iteration in range(n_iterations):
            log.info(f"Contamination correction iteration {iteration + 1} of {n_iterations}")

            per_slit_simuls = list(matched_flat_simuls)

            # Sort fittable slits by decreasing brightness so brighter sources are fitted
            # first.  Their corrected simulations are immediately folded into simul_data
            # so subsequent fainter sources see better contamination within this iteration.
            fittable = [
                k for k in range(len(output_model.slits)) if matched_flat_simuls[k] is not None
            ]
            sort_order = sorted(fittable, key=lambda k: -flat_matched_sum[k])

            # Build simulation from the previous iteration's fitted shapes (flat spectrum for
            # the first iteration).  Update incrementally after each successful fit so that
            # fainter sources benefit from the best available contamination estimate.
            simul_data = _build_simulated_image_from_slits(
                obs.simulated_slits, obs.simulated_image.shape
            )
            success = 0
            for i in sort_order:
                slit = output_model.slits[i]
                matched_flat = matched_flat_simuls[i]

                # Compute the latest contamination estimate
                # For sources brighter than the ith one in the sort order, this will
                # include the polyfit from this order.
                # For fainter sources, it will be whatever was in the previous iteration,
                # i.e., flat-spectrum for the first iteration.
                simul_all_cut = _cut_frame_to_match_slit(simul_data, slit)
                slit.data = original_data[i] - (simul_all_cut - matched_flat.data)

                if _fit_spectral_shape(
                    slit,
                    matched_flat,
                    obs.simulated_slits.slits[good_idxs[i]],
                    polyfit_degree,
                    l2_alpha=l2_alpha,
                    rejection_threshold=rejection_threshold,
                ):
                    success += 1
                    # Immediately rebuild so subsequent (fainter) slits see updated fit
                    simul_data = _build_simulated_image_from_slits(
                        obs.simulated_slits, obs.simulated_image.shape
                    )

            log.info(
                f"Spectral fitting successful for {success} out of {len(output_model.slits)} slits "
                f"in iteration {iteration + 1}. Turn on debug logging for details of failures."
            )

            # Compute per-slit contamination and update corrected data for the next iteration.
            # Always subtract from the original input so errors do not accumulate across iterations.
            contam_cuts = _build_contam(output_model, per_slit_simuls, simul_data, original_data)

            if success == 0:
                log.warning(
                    f"No successful spectral fits in iteration {iteration + 1}. "
                    "Will not continue iterating. Ensure that the background is well subtracted, "
                    "and consider reducing polyfit_degree or increasing l2_alpha to improve fit "
                    "stability."
                )
                break

    # Build output contam_model and simul_slits from the final iteration's results.
    log.info("Creating contamination image for each individual source")
    for i in range(len(per_slit_simuls)):
        this_obs = output_model.slits[i]
        this_simul = per_slit_simuls[i]
        contam_cut = contam_cuts[i]
        if this_simul is not None:
            simul_slit_out = datamodels.SlitModel()
            simul_slit_out.data = this_simul.data
            simul_slit_out.update(this_obs, only="SCI")
            simul_slits.slits.append(simul_slit_out)

        contam_slit = datamodels.SlitModel()
        contam_slit.update(this_obs, only="SCI")
        contam_slit.data = contam_cut
        contam_model.slits.append(contam_slit)

    simul_model = datamodels.ImageModel(data=simul_data)
    simul_model.update(input_model, only="PRIMARY")

    output_model.update(input_model, only="PRIMARY")
    output_model.meta.cal_step.wfss_contam = "COMPLETE"

    return output_model, simul_model, contam_model, simul_slits
