import logging
import warnings

import gwcs
import numpy as np
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.utils.exceptions import AstropyUserWarning
from photutils.background import Background2D, MedianBackground
from scipy.optimize import curve_fit
from stdatamodels.jwst.datamodels import dqflags

from jwst import datamodels
from jwst.assign_wcs import nirspec, AssignWcsStep
from jwst.flatfield import FlatFieldStep
from jwst.clean_flicker_noise.lib import NSClean, NSCleanSubarray
from jwst.lib.basic_utils import LoggingContext
from jwst.msaflagopen import MSAFlagOpenStep
from jwst.ramp_fitting import RampFitStep


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Fixed slit region to mask, for NIRSpec MOS and IFU data
# Values are y start and stop indices, for the edges of the
# region to mask.
NRS_FS_REGION = [922, 1116]


def make_rate(input_model, return_cube=False):
    """
    Make a rate model from a ramp model.

    Parameters
    ----------
    input_model : `~jwst.datamodel.RampModel`
        Input ramp model.

    return_cube : bool, optional
        If set, a CubeModel will be returned, with a separate
        rate for each integration.  Otherwise, an ImageModel is returned
        with the combined rate for the full observation.

    Returns
    -------
    rate_model : `~jwst.datamodel.ImageModel` or `~jwst.datamodel.CubeModel`
        The rate or rateints model.
    """
    # Call the ramp fit step on a copy of the input
    # Use software default values for parameters

    log.info("Creating draft rate file for scene masking")
    step = RampFitStep()
    with LoggingContext(step.log, level=logging.WARNING):
        # Note: the copy is currently needed because ramp fit
        # closes the input model when it's done, and we need
        # it to stay open.
        rate, rateints = step.run(input_model.copy())

    if return_cube:
        output_model = rateints
        rate.close()
        del rate
    else:
        output_model = rate
        rateints.close()
        del rateints

    return output_model


def post_process_rate(input_model, assign_wcs=False, msaflagopen=False, flat_dq=False):
    """
    Post-process the input rate model, as needed.

    Parameters
    ----------
    input_model : `~jwst.datamodel.ImageModel` or `~jwst.datamodel.CubeModel`
        Input rate model.

    assign_wcs : bool, optional
        If set and the input does not already have a WCS assigned,
        the assign_wcs step will be called on the rate model.

    msaflagopen : bool, optional
        If set, the msaflagopen step will be called on the rate model.
        If a WCS is not already present, assign_wcs will be called first.

    flat_dq : bool, optional
        If set, the flat_field step will be run on the input model. DQ
        flags are retrieved from the output and added to the input
        model's DQ array. The rate data is not modified.

    Returns
    -------
    output_model : `~jwst.datamodel.ImageModel` or `~jwst.datamodel.CubeModel`
        The updated model.
    """
    output_model = input_model

    # If needed, assign a WCS
    if (assign_wcs or msaflagopen) and not hasattr(output_model.meta, 'wcs'):
        log.info("Assigning a WCS for scene masking")
        step = AssignWcsStep()
        with LoggingContext(step.log, level=logging.WARNING):
            output_model = step.run(output_model)

    # If needed, flag open MSA shutters
    if msaflagopen:
        log.info("Flagging failed-open MSA shutters for scene masking")
        step = MSAFlagOpenStep()
        with LoggingContext(step.log, level=logging.WARNING):
            output_model = step.run(output_model)

    # If needed, draft a flat correction to retrieve non-science areas
    if flat_dq:
        log.info("Retrieving flat DQ values for scene masking")
        step = FlatFieldStep()
        with LoggingContext(step.log, level=logging.WARNING):
            flat_corrected_model = step.run(output_model)

        # Copy out the flat DQ plane, leave the data as is
        output_model.dq |= flat_corrected_model.dq
        flat_corrected_model.close()
        del flat_corrected_model

    return output_model


def mask_ifu_slices(input_model, mask):
    """
    Flag pixels within IFU slices.

    Find pixels located within IFU slices, according to the WCS,
    and flag them in the mask, so that they do not get used.

    Parameters
    ----------
    input_model : `~jwst.datamodel.JwstDataModel`
        Science data model

    mask : array-like of bool
        2D input mask that will be updated. True indicates background
        pixels to be used. IFU slice regions will be set to False.

    Returns
    -------
    mask : array-like of bool
        2D output mask with additional flags for science pixels
    """
    log.info("Finding slice pixels for an IFU image")

    # Initialize global DQ map to all zero (OK to use)
    dqmap = np.zeros_like(input_model.dq)

    # Get the wcs objects for all IFU slices
    list_of_wcs = nirspec.nrs_ifu_wcs(input_model)

    # Loop over the IFU slices, finding the valid region for each
    for (k, ifu_wcs) in enumerate(list_of_wcs):

        # Construct array indexes for pixels in this slice
        x, y = gwcs.wcstools.grid_from_bounding_box(
            ifu_wcs.bounding_box, step=(1, 1), center=True)

        # Get the world coords for all pixels in this slice;
        # all we actually need are wavelengths
        coords = ifu_wcs(x, y)
        dq = dqmap[..., y.astype(int), x.astype(int)]
        wl = coords[2]

        # set non-NaN wavelength locations as do not use (one)
        valid = ~np.isnan(wl)
        dq = dq[..., valid]
        x = x[..., valid]
        y = y[..., valid]
        dq[:] = 1

        # Copy DQ for this slice into global DQ map
        dqmap[..., y.astype(int), x.astype(int)] = dq

    # Now set all non-zero locations in the mask to False (do not use)
    mask[dqmap == 1] = False

    return mask


def mask_slits(input_model, mask):
    """
    Flag pixels within science regions.

    Find pixels located within MOS or fixed slit footprints
    and flag them in the mask, so that they do not get used.

    Parameters
    ----------
    input_model : `~jwst.datamodel.JwstDataModel`
        Science data model.

    mask : array-like of bool
        2D input mask that will be updated. True indicates background
        pixels to be used. Slit regions will be set to False.

    Returns
    -------
    mask : array-like of bool
        2D output mask with additional flags for slit pixels
    """

    from jwst.extract_2d.nirspec import offset_wcs

    log.info("Finding slit/slitlet pixels")

    # Get the slit-to-msa frame transform from the WCS object
    slit2msa = input_model.meta.wcs.get_transform('slit_frame', 'msa_frame')

    # Loop over the slits, marking all the pixels within each bounding
    # box as False (do not use) in the mask.
    # Note that for 3D masks (TSO mode), all planes will be set to the same value.
    for slit in slit2msa.slits:
        slit_wcs = nirspec.nrs_wcs_set_input(input_model, slit.name)
        xlo, xhi, ylo, yhi = offset_wcs(slit_wcs)
        mask[..., ylo:yhi, xlo:xhi] = False

    return mask


def clip_to_background(image, mask, sigma_lower=3.0, sigma_upper=2.0,
                       fit_histogram=False, lower_half_only=False,
                       verbose=False):
    """
    Flag signal and bad pixels in the image mask.

    Given an image, estimate the background level and a sigma value for the
    mean background.

    The center and sigma may be calculated from a simple sigma-clipped
    median and standard deviation, or may be refined by fitting a Gaussian
    to a histogram of the values.  In that case, the center is the
    Gaussian mean and the sigma is the Gaussian width.

    Pixels above the center + sigma_upper * sigma are assumed to
    have signal; those below this level are assumed to be background pixels.

    Pixels less than center - sigma_lower * sigma are also excluded as bad values.

    The input mask is updated in place.

    Parameters
    ----------
    image : array-like of float
        2D image containing signal and background values.

    mask : array-like of bool
        2D input mask to be updated. True indicates background
        pixels to be used. Regions containing signal or outliers
        will be set to False.

    sigma_lower : float, optional
        The number of standard deviations to use as the lower bound
        for the clipping limit. Values below this limit are marked
        False in the mask.

    sigma_upper : float, optional
        The number of standard deviations to use as the upper bound
        for the clipping limit. Values above this limit are marked
        False in the mask.

    fit_histogram :  bool, optional
        If set, the center value and standard deviation used with
        `sigma_lower` and `sigma_upper` for clipping outliers is derived
        from a Gaussian fit to a histogram of values. Otherwise, the
        center and standard deviation are derived from a simple iterative
        sigma clipping.

    lower_half_only : bool, optional
        If set, the data used to compute the center and standard deviation
        for clipping is the lower half of the distribution only. Values
        below the median are mirrored around the median value to simulate
        a symmetric distribution.  This is intended to account for
        asymmetrical value distributions, with long tails in the upper
        half of the distribution, due to diffuse emission, for example.

    verbose : bool, optional
        If set, DEBUG level messages are issued with details on the
        computed statistics.
    """
    # Sigma limit for basic stats
    sigma_limit = 3.0

    # Check mask for any valid data before proceeding
    if not np.any(mask):
        return

    # Initial iterative sigma clip
    with warnings.catch_warnings():
        warnings.filterwarnings(
            action="ignore", category=AstropyUserWarning)
        mean, median, sigma = sigma_clipped_stats(
            image, mask=~mask, sigma=sigma_limit)
    if fit_histogram:
        center = mean
    else:
        center = median
    if verbose:
        log.debug('From initial sigma clip:')
        log.debug(f'    center: {center:.5g}')
        log.debug(f'    sigma: {sigma:.5g}')

    # If desired, use only the lower half of the data distribution
    if lower_half_only:
        lower_half_idx = mask & (image <= center)
        data_for_stats = np.concatenate(
            ((image[lower_half_idx] - center),
             (center - image[lower_half_idx]))) + center

        # Redo stats on lower half of distribution
        with warnings.catch_warnings():
            warnings.filterwarnings(
                action="ignore", category=AstropyUserWarning)
            mean, median, sigma = sigma_clipped_stats(
                data_for_stats, sigma=sigma_limit)
        if fit_histogram:
            center = mean
        else:
            center = median
        if verbose:
            log.debug('From lower half distribution:')
            log.debug(f'    center: {center:.5g}')
            log.debug(f'    sigma: {sigma:.5g}')
    else:
        data_for_stats = image[mask]

    # Refine sigma and center from a fit to a histogram, if desired
    if fit_histogram:
        try:
            hist, edges = np.histogram(
                data_for_stats, bins=2000,
                range=(center - 4. * sigma, center + 4. * sigma))
        except ValueError:
            log.error('Histogram failed; using clip center and sigma.')
            hist, edges, param_opt = None, None, None

        if hist is not None:
            values = (edges[1:] + edges[0:-1]) / 2.
            ind = np.argmax(hist)
            mode_estimate = values[ind]

            # Fit a Gaussian profile to the histogram
            def gaussian(x, g_amp, g_mean, g_sigma):
                return g_amp * np.exp(-0.5 * ((x - g_mean) / g_sigma) ** 2)

            param_start = (hist[ind], mode_estimate, sigma)
            bounds = [(0, values[0], 0),
                      (np.inf, values[-1], values[-1] - values[0])]
            try:
                param_opt, _ = curve_fit(gaussian, values, hist, p0=param_start,
                                         bounds=bounds)
            except RuntimeError:
                log.error('Gaussian fit failed; using clip center and sigma.')
                param_opt = None

            if verbose:
                log.debug('From histogram:')
                log.debug(f'    mode estimate: {mode_estimate:.5g}')
                log.debug(f'    range of values in histogram: '
                          f'{values[0]:.5g} to {values[-1]:.5g}')

        if verbose:
            log.debug('Gaussian fit results:')
        if param_opt is None:
            if verbose:
                log.debug('    (fit failed)')
        else:
            if verbose:
                log.debug(f'    peak: {param_opt[0]:.5g}')
                log.debug(f'    center: {param_opt[1]:.5g}')
                log.debug(f'    sigma: {param_opt[2]:.5g}')
            center = param_opt[1]
            sigma = param_opt[2]

    # Set limits from center and sigma
    background_lower_limit = center - sigma_lower * sigma
    background_upper_limit = center + sigma_upper * sigma
    if verbose:
        log.debug(f'Mask limits: {background_lower_limit:.5g} '
                  f'to {background_upper_limit:.5g}')

    # Clip bad values
    bad_values = image < background_lower_limit
    mask[bad_values] = False

    # Clip signal (> N sigma)
    signal = image > background_upper_limit
    mask[signal] = False


def create_mask(input_model, mask_science_regions=False,
                n_sigma=2.0, fit_histogram=False, single_mask=False):
    """
    Create a mask identifying background pixels.

    Parameters
    ----------
    input_model : `~jwst.datamodel.JwstDataModel`
        Science data model, containing rate data with all necessary
        pre-processing already performed.

    mask_science_regions : bool, optional
        For NIRSpec, mask regions of the image defined by WCS bounding
        boxes for slits/slices, as well as any regions known to be
        affected by failed-open MSA shutters. This requires the
        `assign_wcs` and `msaflagopen` steps to have been run on the
        input_model.
        For MIRI imaging, mask regions of the detector not used for science.
        This requires that NON_SCIENCE flags are set in the DQ array
        for the input_model.

    n_sigma : float, optional
        Sigma threshold for masking outliers.

    fit_histogram : bool, optional
        If set, the 'sigma' used with `n_sigma` for clipping outliers
        is derived from a Gaussian fit to a histogram of values.
        Otherwise, a simple iterative sigma clipping is performed.

    single_mask : bool, optional
        If set, a single mask will be created, regardless of
        the number of input integrations. Otherwise, the mask will
        be a 3D cube, with one plane for each integration.

    Returns
    -------
    mask : array-like of bool
        2D or 3D image mask
    """
    exptype = input_model.meta.exposure.type.lower()

    # Initialize mask to all True. Subsequent operations will mask
    # out pixels that contain signal.
    mask = np.full(input_model.dq.shape, True)

    # If IFU, mask all pixels contained in the IFU slices
    if exptype == 'nrs_ifu' and mask_science_regions:
        mask = mask_ifu_slices(input_model, mask)

    # If MOS or FS, mask all pixels affected by open slitlets
    if (exptype in ['nrs_fixedslit', 'nrs_brightobj', 'nrs_msaspec']
            and mask_science_regions):
        mask = mask_slits(input_model, mask)

    # If IFU or MOS, mask pixels affected by failed-open shutters
    if mask_science_regions and exptype in ['nrs_ifu', 'nrs_msaspec']:
        open_pix = input_model.dq & dqflags.pixel['MSA_FAILED_OPEN']
        mask[open_pix > 0] = False

    # If MIRI imaging, mask the non-science regions:
    # they contain irrelevant emission
    if mask_science_regions and exptype in ['mir_image']:
        non_science = (input_model.dq & dqflags.pixel['NON_SCIENCE']) > 1
        mask[non_science] = False

    # Mask any NaN pixels or exactly zero value pixels
    no_data_pix = np.isnan(input_model.data) | (input_model.data == 0)
    mask[no_data_pix] = False

    # If IFU or MOS, mask the fixed-slit area of the image; uses hardwired indexes
    if mask_science_regions:
        if exptype == 'nrs_ifu':
            log.info("Masking the fixed slit region for IFU data.")
            mask[..., NRS_FS_REGION[0]:NRS_FS_REGION[1], :] = False
        elif exptype == 'nrs_msaspec':
            # check for any slits defined in the fixed slit quadrant:
            # if there is nothing there of interest, mask the whole FS region
            try:
                slit2msa = input_model.meta.wcs.get_transform('slit_frame', 'msa_frame')
                is_fs = [s.quadrant == 5 for s in slit2msa.slits]
            except (AttributeError, ValueError, TypeError):
                log.warning('Slit to MSA transform not found.')
                is_fs = [False]
            if not any(is_fs):
                log.info("Masking the fixed slit region for MOS data.")
                mask[..., NRS_FS_REGION[0]:NRS_FS_REGION[1], :] = False
            else:
                log.info("Fixed slits found in MSA definition; "
                         "not masking the fixed slit region for MOS data.")

    # Mask outliers and signal using sigma clipping stats.
    # For 3D data, loop over each integration separately.
    if input_model.data.ndim == 3:
        for i in range(input_model.data.shape[0]):
            clip_to_background(
                input_model.data[i], mask[i],
                sigma_upper=n_sigma, fit_histogram=fit_histogram, verbose=True)
    else:
        clip_to_background(
            input_model.data, mask,
            sigma_upper=n_sigma, fit_histogram=fit_histogram, verbose=True)

    # Reduce the mask to a single plane if needed
    if single_mask and mask.ndim == 3:
        mask = np.all(mask, axis=0)

    # Return the mask
    return mask


def background_level(image, mask, background_method='median',
                     background_box_size=None):
    """
    Fit a low-resolution background level.

    Parameters
    ----------
    image : array-like of float
        The 2D image containing the background to fit.

    mask : array-like of bool
        The mask that indicates which pixels are to be used in fitting.
        True indicates a background pixel.

    background_method : {'median', 'model', None}, optional
        If 'median', the preliminary background to remove and restore
        is a simple median of the background data.  If 'model', the
        background data is fit with a low-resolution model via
        `~photutils.background.Background2D`.  If None, the background
        value is 0.0.

    background_box_size : tuple of int, optional
        Box size for the data grid used by `Background2D` when
        `background_method` = 'model'. For best results, use a box size
        that evenly divides the input image shape. Defaults to 32x32 if
        not provided.

    Returns
    -------
    background : float or array-like of float
        The background level: a single value, if `background_method`
        is 'median' or None, or an array matching the input image size
        if `background_method` is 'model'.
    """
    if background_method is None:
        background = 0.0

    else:
        # Sigma limit for basic stats
        sigma_limit = 3.0

        # Flag more signal in the background subtracted image,
        # with sigma set by the lower half of the distribution only
        clip_to_background(
            image, mask, sigma_lower=sigma_limit, sigma_upper=sigma_limit,
            lower_half_only=True)

        if background_method == 'model':
            sigma_clip_for_bkg = SigmaClip(sigma=sigma_limit, maxiters=5)
            bkg_estimator = MedianBackground()

            if background_box_size is None:
                background_box_size = [32, 32]
            box_division_remainder = (image.shape[0] % background_box_size[0],
                                      image.shape[1] % background_box_size[1])
            if not np.allclose(box_division_remainder, 0):
                log.warning(f'Background box size {background_box_size} '
                            f'does not divide evenly into the image '
                            f'shape {image.shape}.')

            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        action="ignore", category=AstropyUserWarning)
                    bkg = Background2D(
                        image, box_size=background_box_size, filter_size=(5, 5),
                        mask=~mask, sigma_clip=sigma_clip_for_bkg,
                        bkg_estimator=bkg_estimator)
                background = bkg.background
            except ValueError:
                log.error('Background fit failed, using median value.')
                background = np.nanmedian(image[mask])
        else:
            background = np.nanmedian(image[mask])
    return background


def fft_clean_full_frame(image, mask, detector):
    """
    Fit and remove background noise in frequency space for a full-frame image.

    Parameters
    ----------
    image : array-like of float
        The image to be cleaned.

    mask : array-like of bool
        The mask that indicates which pixels are to be used in fitting.

    detector : str
        The name of the detector from which the data originate.

    Returns
    -------
    cleaned_image : array-like of float
        The cleaned image.
    """

    # Instantiate the cleaner
    cleaner = NSClean(detector, mask)

    # Clean the image
    try:
        cleaned_image = cleaner.clean(image, buff=True)
    except np.linalg.LinAlgError:
        log.warning("Error cleaning image; step will be skipped")
        return None

    return cleaned_image


def fft_clean_subarray(image, mask, detector):
    """
    Fit and remove background noise in frequency space for a subarray image.

    Parameters
    ----------
    image : array-like of float
        The image to be cleaned.

    mask : array-like of bool
        The mask that indicates which pixels are to be used in fitting.

    detector : str
        The name of the detector from which the data originate.

    Returns
    -------
    cleaned_image : array-like of float
        The cleaned image.
    """

    # Flip the image to detector coords. NRS1 requires a transpose
    # of the axes, while NRS2 requires a transpose and flip.
    if detector == "NRS1":
        image = image.transpose()
        mask = mask.transpose()
    else:
        image = image.transpose()[::-1]
        mask = mask.transpose()[::-1]

    # Instantiate the cleaner
    cleaner = NSCleanSubarray(image, mask)

    # Clean the image
    try:
        cleaned_image = cleaner.clean()
    except np.linalg.LinAlgError:
        log.warning("Error cleaning image; step will be skipped")
        return None

    # Restore the cleaned image to the science frame
    if detector == "NRS1":
        cleaned_image = cleaned_image.transpose()
    else:
        cleaned_image = cleaned_image[::-1].transpose()

    return cleaned_image


def median_clean(image, mask, axis_to_correct):
    """
    Fit and remove background noise via median values along one image axis.

    Parameters
    ----------
    image : array-like of float
        The image to be cleaned.

    mask : array-like of bool
        The mask that indicates which pixels are to be used in fitting.
        True indicates a background pixel.

    axis_to_correct : int
        For NIR detectors, the axis to correct should be the detector slow
        readout direction.  Values expected are 1 or 2, following the JWST
        datamodel definition (meta.subarray.slowaxis).  For MIRI, flicker
        noise appears along the vertical direction, so `axis_to_correct`
        should be set to 1 (median along the y-axis).

    Returns
    -------
    cleaned_image : array-like of float
        The cleaned image.
    """
    # Masked median along slow axis
    masked_image = np.ma.array(image, mask=~mask)
    stripes = np.ma.median(masked_image, axis=(axis_to_correct - 1), keepdims=True)
    stripes = np.ma.filled(stripes, fill_value=0.0)

    # Remove median stripes
    corrected_image = image - stripes

    return corrected_image


def do_correction(input_model, fit_method='median', background_method='median',
                  background_box_size=None, background_from_rate=False,
                  mask_science_regions=False, n_sigma=2.0, fit_histogram=False,
                  single_mask=True, user_mask=None, save_mask=False,
                  save_background=False, save_noise=False):
    """
    Apply the 1/f noise correction.

    Parameters
    ----------
    input_model : `~jwst.datamodel.JwstDataModel`
        Science data to be corrected.

    fit_method : {'fft', 'median'}, optional
        The algorithm to use to fit background noise.

    background_method : {'median', 'model', None}, optional
        If 'median', the preliminary background to remove and restore
        is a simple median of the background data.  If 'model', the
        background data is fit with a low-resolution model via
        `~photutils.background.Background2D`.  If None, the background
        value is 0.0.

    background_box_size : tuple of int, optional
        Box size for the data grid used by `Background2D` when
        `background_method` = 'model'. For best results, use a box size
        that evenly divides the input image shape.

    background_from_rate : bool, optional
        If set, and the input data is a ramp model, the background will be
        fit to the rate images instead of the individual group differences.
        The preliminary background subtracted from each diff before fitting
        noise is then rate background * group time.

    mask_science_regions : bool, optional
        Mask slit/slice regions defined in WCS. Implemented only for
        NIRSpec science data modes.

    n_sigma : float, optional
        N-sigma rejection level for finding outliers.

    fit_histogram : bool, optional
        If set, the 'sigma' used with `n_sigma` for clipping outliers
        is derived from a Gaussian fit to a histogram of values.
        Otherwise, a simple iterative sigma clipping is performed.

    single_mask : bool, optional
        If set, a single mask will be created, regardless of
        the number of input integrations. Otherwise, the mask will
        be a 3D cube, with one plane for each integration.

    user_mask : str or None, optional
        Path to user-supplied mask image.

    save_mask : bool, optional
        Switch to indicate whether the mask should be saved.

    save_background : bool, optional
        Switch to indicate whether the fit background should be saved.

    save_noise : bool, optional
        Switch to indicate whether the fit noise should be saved.

    Returns
    -------
    output_model : `~jwst.datamodel.JwstDataModel`
        Corrected data.

    mask_model : `~jwst.datamodel.JwstDataModel`
        Pixel mask to be saved or None.

    background_model : `~jwst.datamodel.JwstDataModel`
        Background model to be saved or None.

    noise_model : `~jwst.datamodel.JwstDataModel`
        Background model to be saved or None.

    status : {'COMPLETE', 'SKIPPED'}
        Completion status.  If errors were encountered, status = 'SKIPPED'
        and the output data matches the input data.  Otherwise,
        status = 'COMPLETE'.
    """
    # Track the completion status, for various failure conditions
    status = 'SKIPPED'

    detector = input_model.meta.instrument.detector.upper()
    subarray = input_model.meta.subarray.name.upper()
    exp_type = input_model.meta.exposure.type
    log.info(f'Input exposure type is {exp_type}, detector={detector}')

    # Check for a valid input that we can work on
    message = None
    miri_allowed = ['MIR_IMAGE']
    if exp_type.startswith('MIR') and exp_type not in miri_allowed:
        message = f"EXP_TYPE {exp_type} is not supported"

    nsclean_allowed = ['NRS_MSASPEC', 'NRS_IFU', 'NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']
    if fit_method == 'fft':
        if exp_type not in nsclean_allowed:
            message = f"Fit method 'fft' cannot be applied to exp_type {exp_type}"
        elif subarray == 'ALLSLITS':
            message = f"Fit method 'fft' cannot be applied to subarray {subarray}"

    if message is not None:
        log.warning(message)
        log.warning("Step will be skipped")
        return input_model, None, None, None, status

    output_model = input_model.copy()

    # Get axis to correct, by instrument
    if exp_type.startswith('MIR'):
        # MIRI doesn't have f-noise, but it does have a vertical flickering.
        # Set the axis for median correction to the y-axis.
        axis_to_correct = 1
    else:
        # For NIR detectors, the axis for median correction is the
        # detector slowaxis.
        axis_to_correct = abs(input_model.meta.subarray.slowaxis)

    # Standardize background arguments
    if str(background_method).lower() == 'none':
        background_method = None
    if background_method is None:
        background_from_rate = False

    # Make a rate file if needed
    if user_mask is None or background_from_rate:
        if isinstance(input_model, datamodels.RampModel):
            image_model = make_rate(input_model, return_cube=(not single_mask))
        else:
            # input is already a rate file
            image_model = input_model

        # If needed, assign a WCS to the rate file,
        # flag open MSA shutters, or retrieve flat DQ flags
        assign_wcs = exp_type.startswith('NRS')
        flag_open = (exp_type in ['NRS_IFU', 'NRS_MSASPEC'])
        flat_dq = (exp_type in ['MIR_IMAGE'])
        if mask_science_regions:
            image_model = post_process_rate(
                image_model, assign_wcs=assign_wcs,
                msaflagopen=flag_open, flat_dq=flat_dq)
    else:
        image_model = input_model

    # Check for a user-supplied mask image. If provided, use it.
    if user_mask is not None:
        mask_model = datamodels.open(user_mask)
        background_mask = (mask_model.data.copy()).astype(np.bool_)
    else:
        # Create the pixel mask that will be used to indicate which pixels
        # to include in the 1/f noise measurements. Basically, we're setting
        # all illuminated pixels to False, so that they do not get used, and
        # setting all unilluminated pixels to True (so they DO get used).
        log.info("Creating mask")
        background_mask = create_mask(
            image_model,
            mask_science_regions=mask_science_regions,
            n_sigma=n_sigma, fit_histogram=fit_histogram,
            single_mask=single_mask)

        # Store the mask image in a model, if requested
        if save_mask:
            if len(background_mask.shape) == 3:
                mask_model = datamodels.CubeModel(data=background_mask)
            else:
                mask_model = datamodels.ImageModel(data=background_mask)

            # Copy metadata from input
            mask_model.update(output_model)
        else:
            mask_model = None

    log.info(f"Cleaning image {input_model.meta.filename}")

    # Setup for handling 2D, 3D, or 4D inputs
    image_shape = input_model.data.shape[-2:]
    ndim = input_model.data.ndim
    group_time = 1.0
    if ndim == 2:
        nints = 1
        ngroups = 1
    else:
        nints = input_model.data.shape[0]
        if ndim == 3:
            ngroups = 1
        else:
            ngroups = input_model.data.shape[1] - 1
            group_time = input_model.meta.exposure.group_time

        # Check for 3D mask
        if background_mask.ndim == 2:
            if nints > 1:
                log.info("Data has multiple integrations, but mask is 2D: "
                         "the same mask will be used for all integrations.")
        elif background_mask.shape[0] != nints:
            log.warning("Mask does not match data shape. Step will be skipped.")
            return output_model, None, None, None, status

    # Check that mask matches 2D data shape
    if background_mask.shape[-2:] != image_shape:
        log.warning("Mask does not match data shape. Step will be skipped.")
        return output_model, None, None, None, status

    # Get the background level from the draft rate file if desired
    if ndim > 3 and background_from_rate:
        if single_mask:
            # One rate image, use for all integrations
            rate_background = background_level(
                image_model.data, background_mask,
                background_method=background_method,
                background_box_size=background_box_size)
            rate_background *= group_time
        else:
            # Compute background for each integration separately
            rate_background = []
            for integ in range(nints):
                rbg = background_level(
                    image_model.data[integ], background_mask[integ],
                    background_method=background_method,
                    background_box_size=background_box_size)
                rate_background.append(rbg * group_time)
    else:
        background_from_rate = False
        rate_background = None

    # Close the draft rate model if created - it is no longer needed.
    if image_model is not input_model:
        image_model.close()
        del image_model

    # Make a background cube for saving, if desired
    if save_background:
        background_to_save = np.zeros_like(input_model.data)
    else:
        background_to_save = None

    # Loop over integrations and groups (even if there's only 1)
    for i in range(nints):
        log.debug(f"Working on integration {i + 1}")
        for j in range(ngroups):
            log.debug(f"Working on group {j + 1}")

            # Copy the scene mask, for further flagging
            if background_mask.ndim == 3:
                mask = background_mask[i].copy()
            else:
                mask = background_mask.copy()

            # Get the relevant image data
            if ndim == 2:
                image = input_model.data
            elif ndim == 3:
                image = input_model.data[i]
            else:
                # Ramp data input:
                # subtract the current group from the next one
                image = input_model.data[i, j+1] - input_model.data[i, j]
                dq = input_model.groupdq[i, j+1]

                # For ramp data, mask any DNU and JUMP pixels
                dnu = (dq & dqflags.group['DO_NOT_USE']) > 0
                mask[dnu] = False
                jump = (dq & dqflags.group['JUMP_DET']) > 0
                mask[jump] = False

            # Make sure data is float32
            image = image.astype(np.float32)

            # Find and replace/mask NaNs
            nan_pix = np.isnan(image)
            image[nan_pix] = 0.0
            mask[nan_pix] = False

            # If there's no good data remaining in the mask,
            # skip this image
            if not np.any(mask):
                log.warning(f'No usable data in integration {i+1}, group {j+1}. '
                            f'Skipping correction for this image.')
                continue

            # Fit and remove a background level
            if str(background_method).lower() == 'none':
                background = 0.0
                bkg_sub = image
            else:
                if background_from_rate:
                    if single_mask:
                        background = rate_background
                    else:
                        background = rate_background[i]
                else:
                    background = background_level(
                        image, mask, background_method=background_method,
                        background_box_size=background_box_size)
                log.debug(f'Background level: {np.nanmedian(background):.5g}')
                bkg_sub = image - background

                # Flag more signal in the background subtracted image,
                # with sigma set by the lower half of the distribution only
                clip_to_background(
                    bkg_sub, mask, sigma_lower=n_sigma,
                    sigma_upper=n_sigma, lower_half_only=True)

            if fit_method == 'fft':
                if bkg_sub.shape == (2048, 2048):
                    # Clean a full-frame image
                    cleaned_image = fft_clean_full_frame(bkg_sub, mask, detector)
                else:
                    # Clean a subarray image
                    cleaned_image = fft_clean_subarray(bkg_sub, mask, detector)
            else:
                cleaned_image = median_clean(bkg_sub, mask, axis_to_correct)

            # Check for failure
            if cleaned_image is None:
                log.error(f'Cleaning failed for integration {i + 1}, group {j + 1}')

                # Restore input data to make sure any partial changes are thrown away
                output_model.data = input_model.data.copy()
                return output_model, None, None, None, status
            else:
                # Restore the background level
                cleaned_image += background

                # Restore NaNs in cleaned image
                cleaned_image[nan_pix] = np.nan

                # Store the cleaned image in the output model
                if ndim == 2:
                    output_model.data = cleaned_image
                    if save_background:
                        background_to_save[:] = background
                elif ndim == 3:
                    output_model.data[i] = cleaned_image
                    if save_background:
                        background_to_save[i] = background
                else:
                    # add the cleaned data diff to the previously cleaned group,
                    # rather than the noisy input group
                    output_model.data[i, j+1] = output_model.data[i, j] + cleaned_image
                    if save_background:
                        background_to_save[i, j+1] = background

    # Store the background image in a model, if requested
    if save_background:
        if background_to_save.ndim == 4:
            background_model = datamodels.RampModel(data=background_to_save)
        elif background_to_save.ndim == 3:
            background_model = datamodels.CubeModel(data=background_to_save)
        else:
            background_model = datamodels.ImageModel(data=background_to_save)

        # Copy metadata from input
        background_model.update(output_model)
    else:
        background_model = None

    # For the noise image, diff the input and output models
    if save_noise:
        noise_data = output_model.data - input_model.data
        if input_model.data.ndim == 4:
            noise_model = datamodels.RampModel(data=noise_data)
        elif input_model.data.ndim == 3:
            noise_model = datamodels.CubeModel(data=noise_data)
        else:
            noise_model = datamodels.ImageModel(data=noise_data)

        # Copy metadata from input
        noise_model.update(output_model)
    else:
        noise_model = None

    # Set completion status
    status = 'COMPLETE'

    return output_model, mask_model, background_model, noise_model, status
