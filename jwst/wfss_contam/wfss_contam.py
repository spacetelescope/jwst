"""Top-level module for WFSS contamination correction."""

import logging
import multiprocessing

import numpy as np
from astropy.table import Table
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.transforms.models import (
    NIRCAMBackwardGrismDispersion,
    NIRISSBackwardGrismDispersion,
)

from jwst.wfss_contam.observations import Observation
from jwst.wfss_contam.sens1d import get_photom_data

log = logging.getLogger(__name__)

__all__ = ["contam_corr"]


def determine_multiprocessing_ncores(max_cores, num_cores):
    """
    Determine the number of cores to use for multiprocessing.

    Parameters
    ----------
    max_cores : str or int
        Number of cores to use for multiprocessing. If set to 'none'
        (the default), then no multiprocessing will be done. The other
        allowable string values are 'quarter', 'half', and 'all', which indicate
        the fraction of cores to use for multi-proc. The total number of
        cores includes the SMT cores (Hyper Threading for Intel).
        If an integer is provided, it will be the exact number of cores used.
    num_cores : int
        Number of cores available on the machine

    Returns
    -------
    ncpus : int
        Number of cores to use for multiprocessing
    """
    match max_cores:
        case "none":
            return 1
        case None:
            return 1
        case "quarter":
            return num_cores // 4 or 1
        case "half":
            return num_cores // 2 or 1
        case "all":
            return num_cores
        case int():
            if max_cores <= num_cores and max_cores > 0:
                return max_cores
            log.warning(
                f"Requested {max_cores} cores exceeds the number of cores available "
                "on this machine ({num_cores}). Using all available cores."
            )
            return num_cores
        case _:
            raise ValueError(f"Invalid value for max_cores: {max_cores}")


class UnmatchedSlitIDError(Exception):
    """Exception raised when a slit ID is not found in the list of simulated slits."""

    pass


def _find_matching_simul_slit(slit, simul_slit_sids, simul_slit_orders):
    """
    Find the index of the matching simulated slit in the list of simulated slits.

    Parameters
    ----------
    slit : `~jwst.datamodels.SlitModel`
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
    slit : `~jwst.datamodels.SlitModel`
        Source slit model

    Returns
    -------
    cutout : 2D array
        Contamination image cutout that matches the extent of the source slit
    """
    x1 = slit.xstart
    y1 = slit.ystart
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
    slit0 : `~jwst.datamodels.SlitModel`
        Slit model for the first slit, which is used as reference.
    slit1 : `~jwst.datamodels.SlitModel`
        Slit model for the second slit, which is reshaped to match slit0.

    Returns
    -------
    slit0, slit1 : `~jwst.datamodels.SlitModel`
        Reshaped slit models slit0, slit1.
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

    backplane1[i0:i1, j0:j1] = data1[i0:i1, j0:j1]

    slit1.data = backplane1
    slit1.xstart = slit0.xstart
    slit1.ystart = slit0.ystart
    slit1.xsize = slit0.xsize
    slit1.ysize = slit0.ysize

    return slit0, slit1


def match_backplane_encompass_both(slit0, slit1):
    """
    Put data from the two slits into a common backplane, encompassing both.

    Slits are zero-padded where their new extent does not overlap with the original data.

    Parameters
    ----------
    slit0, slit1 : `~jwst.datamodels.SlitModel`
        Slit model for the first and second slit.

    Returns
    -------
    slit0, slit1 : `~jwst.datamodels.SlitModel`
        Reshaped slit models slit0, slit1.
    """
    data0 = slit0.data
    data1 = slit1.data

    shape = (max(data0.shape[0], data1.shape[0]), max(data0.shape[1], data1.shape[1]))
    xmin = min(slit0.xstart, slit1.xstart)
    ymin = min(slit0.ystart, slit1.ystart)
    shape = (
        max(
            slit0.xsize + slit0.xstart - xmin,
            slit1.xsize + slit1.xstart - xmin,
        ),
        max(
            slit0.ysize + slit0.ystart - ymin,
            slit1.ysize + slit1.ystart - ymin,
        ),
    )
    x0 = slit0.xstart - xmin
    y0 = slit0.ystart - ymin
    x1 = slit1.xstart - xmin
    y1 = slit1.ystart - ymin

    backplane0 = np.zeros(shape).T
    backplane0[y0 : y0 + data0.shape[0], x0 : x0 + data0.shape[1]] = data0
    backplane1 = np.zeros(shape).T
    backplane1[y1 : y1 + data1.shape[0], x1 : x1 + data1.shape[1]] = data1

    slit0.data = backplane0
    slit1.data = backplane1
    for slit in [slit0, slit1]:
        slit.xstart = xmin
        slit.ystart = ymin
        slit.xsize = shape[0]
        slit.ysize = shape[1]

    return slit0, slit1


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


def contam_corr(
    input_model,
    waverange,
    photom,
    max_cores,
    orders=None,
    magnitude_limit=None,
    max_pixels_per_chunk=5e4,
    oversample_factor=2,
):
    """
    Correct contamination in WFSS spectral cutouts.

    Parameters
    ----------
    input_model : `~jwst.datamodels.MultiSlitModel`
        Input data model containing 2D spectral cutouts
    waverange : `~jwst.datamodels.WavelengthrangeModel`
        Wavelength range reference file model
    photom : `~jwst.datamodels.NrcWfssPhotomModel` or `~jwst.datamodels.NisWfssPhotomModel`
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

    Returns
    -------
    output_model : `~jwst.datamodels.MultiSlitModel`
        A copy of the input_model that has been corrected
    simul_model : `~jwst.datamodels.ImageModel`
        Full-frame simulated image of the grism exposure
    contam_model : `~jwst.datamodels.MultiSlitModel`
        Contamination estimate images for each source slit
    """
    num_cores = multiprocessing.cpu_count()
    ncpus = determine_multiprocessing_ncores(max_cores, num_cores)

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
        phot_per_lam = False
    else:
        filter_name = filter_kwd
        phot_per_lam = True

    # Read the source catalog to perform magnitude-based source selection later
    # mag limit will be scaled according to order 1 sensitivity
    if magnitude_limit is not None:
        source_catalog = Table.read(input_model.meta.source_catalog, format="ascii.ecsv")
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
        phot_per_lam=phot_per_lam,
    )

    no_sources = True
    for order in spec_orders:
        # Load lists of wavelength ranges and flux cal info
        wavelength_range = waverange.get_wfss_wavelength_range(filter_name, [order])
        wmin = wavelength_range[order][0]
        wmax = wavelength_range[order][1]
        log.debug(f"wmin={wmin}, wmax={wmax} for order {order}")
        sens_waves, sens_response = get_photom_data(photom, filter_kwd, pupil_kwd, order)

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
        obs.disperse_order(order, wmin, wmax, sens_waves, sens_response, selected_ids)

    if no_sources:
        log.error(
            f"No sources found that met the magnitude limit {magnitude_limit}. Step will be SKIPPED"
        )
        return input_model, None, None, None

    # Initialize the full-frame simulated grism image
    simul_model = datamodels.ImageModel(data=obs.simulated_image)
    simul_model.update(input_model, only="PRIMARY")

    simul_slit_sids = [slit.source_id for slit in obs.simulated_slits.slits]
    simul_slit_orders = [slit.meta.wcsinfo.spectral_order for slit in obs.simulated_slits.slits]

    # Initialize output multislitmodel
    output_model = datamodels.MultiSlitModel()

    # Copy over matching slits
    good_slits = [
        datamodels.SlitModel(slit.instance).copy()
        for slit in input_model.slits
        if slit.source_id in obs.source_ids
    ]
    output_model.slits.extend(good_slits)

    # Loop over all slits/sources to subtract contaminating spectra
    log.info("Creating contamination image for each individual source")
    contam_model = datamodels.MultiSlitModel()
    contam_model.update(input_model, only="PRIMARY")
    simul_slits = datamodels.MultiSlitModel()
    simul_slits.update(input_model, only="PRIMARY")
    for slit in output_model.slits:
        try:
            good_idx = _find_matching_simul_slit(slit, simul_slit_sids, simul_slit_orders)
            this_simul = obs.simulated_slits.slits[good_idx]
            slit, this_simul = match_backplane_prefer_first(slit, this_simul)
            simul_all_cut = _cut_frame_to_match_slit(obs.simulated_image, slit)
            contam_cut = simul_all_cut - this_simul.data
            simul_slits.slits.append(this_simul)

        except (UnmatchedSlitIDError, SlitOverlapError) as e:
            log.warning(e)
            contam_cut = np.zeros_like(slit.data)

        contam_slit = datamodels.SlitModel()
        contam_slit.data = contam_cut
        contam_model.slits.append(contam_slit)

        # Subtract the contamination from the source slit
        slit.data -= contam_cut

    output_model.update(input_model, only="PRIMARY")
    output_model.meta.cal_step.wfss_contam = "COMPLETE"
    seg_model.close()

    return output_model, simul_model, contam_model, simul_slits
