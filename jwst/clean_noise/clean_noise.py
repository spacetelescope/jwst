import logging

import gwcs
import numpy as np
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.background import Background2D, MedianBackground
from scipy.optimize import curve_fit
from stdatamodels.jwst.datamodels import dqflags

from jwst import datamodels
from jwst.assign_wcs import nirspec, AssignWcsStep
from jwst.clean_noise.lib import NSClean, NSCleanSubarray
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
    # Note: the copy is currently needed because ramp fit
    # closes the input model when it's done, and we need
    # it to stay open.
    log.info("Creating draft rate file for scene masking")
    step = RampFitStep()
    with LoggingContext(step.log, level=logging.WARNING):
        # Use software default values for parameters
        rate, rateints = step.run(input_model.copy())

    if return_cube:
        output_model = rateints
        rate.close()
    else:
        output_model = rate
        rateints.close()

    return output_model


def assign_wcs_to_rate(input_model, msaflagopen=False):
    """
    Assign a WCS to the input rate model.

    Parameters
    ----------
    input_model : `~jwst.datamodel.ImageModel` or `~jwst.datamodel.CubeModel`
        Input rate model.

    msaflagopen : bool, optional
        If set, the msaflagopen step will be additionally be called
        on the rate model before returning.

    Returns
    -------
    output_model : `~jwst.datamodel.ImageModel` or `~jwst.datamodel.CubeModel`
        The updated model.
    """
    log.info("Assigning a WCS for scene masking")
    step = AssignWcsStep()
    with LoggingContext(step.log, level=logging.WARNING):
        output_model = step.run(input_model)

    # If needed, flag open MSA shutters
    if msaflagopen:
        log.info("Flagging failed-open MSA shutters for scene masking")
        step = MSAFlagOpenStep()
        with LoggingContext(step.log, level=logging.WARNING):
            output_model = step.run(output_model)

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
        2D input mask that will be updated

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
        x, y = gwcs.wcstools.grid_from_bounding_box(ifu_wcs.bounding_box,
                                                    step=(1, 1),
                                                    center=True)
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
        Science data model

    mask : array-like of bool
        2D input mask that will be updated

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
    image
    mask
    sigma_lower
    sigma_upper
    fit_histogram
    lower_half_only
    verbose
    """
    # Sigma limit for basic stats
    sigma_limit = 3.0

    # Initial iterative sigma clip
    mean, median, sigma = sigma_clipped_stats(image, mask=~mask, sigma=sigma_limit)
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
        lower_half_idx = mask & (image < center)
        data_for_stats = np.concatenate(
            ((image[lower_half_idx] - center),
             (center - image[lower_half_idx]))) + center

        # Redo stats on lower half of distribution
        mean, median, sigma = sigma_clipped_stats(data_for_stats, sigma=sigma_limit)
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
        hist, edges = np.histogram(data_for_stats, bins=2000,
                                   range=(center - 4. * sigma, center + 4. * sigma))
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

    return


def create_mask(input_model, mask_spectral_regions=False,
                n_sigma=2.0, fit_histogram=False, single_mask=False):
    """
    Create a mask identifying background pixels.

    Parameters
    ----------
    input_model : `~jwst.datamodel.JwstDataModel`
        Science data model.

    mask_spectral_regions : bool, optional
        Mask slit/slice regions defined in WCS. Implemented
        only for NIRSpec science modes.

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

    # make a rate file if needed
    flag_open = (exptype in ['nrs_ifu', 'nrs_msaspec'])
    if isinstance(input_model, datamodels.RampModel):
        image_model = make_rate(input_model, return_cube=(not single_mask))

        # if needed, also assign a WCS
        if mask_spectral_regions:
            image_model = assign_wcs_to_rate(image_model, msaflagopen=flag_open)
    else:
        # input is already a rate file
        image_model = input_model

        # if needed, assign a WCS
        if mask_spectral_regions and not hasattr(image_model.meta, 'wcs'):
            image_model = assign_wcs_to_rate(image_model, msaflagopen=flag_open)

    # Initialize mask to all True. Subsequent operations will mask
    # out pixels that contain signal.
    mask = np.full(image_model.dq.shape, True)

    # If IFU, mask all pixels contained in the IFU slices
    if exptype == 'nrs_ifu' and mask_spectral_regions:
        mask = mask_ifu_slices(image_model, mask)

    # If MOS or FS, mask all pixels affected by open slitlets
    if (exptype in ['nrs_fixedslit', 'nrs_brightobj', 'nrs_msaspec']
            and mask_spectral_regions):
        mask = mask_slits(image_model, mask)

    # If IFU or MOS, mask pixels affected by failed-open shutters
    if mask_spectral_regions and exptype in ['nrs_ifu', 'nrs_msaspec']:
        open_pix = image_model.dq & dqflags.pixel['MSA_FAILED_OPEN']
        mask[open_pix > 0] = False

    # Mask any NaN pixels or exactly zero value pixels
    no_data_pix = np.isnan(image_model.data) | (image_model.data == 0)
    mask[no_data_pix] = False

    # If IFU or MOS, mask the fixed-slit area of the image; uses hardwired indexes
    if mask_spectral_regions:
        if exptype == 'nrs_ifu':
            log.info("Masking the fixed slit region for IFU data.")
            mask[..., NRS_FS_REGION[0]:NRS_FS_REGION[1], :] = False
        elif exptype == 'nrs_msaspec':
            # check for any slits defined in the fixed slit quadrant:
            # if there is nothing there of interest, mask the whole FS region
            slit2msa = image_model.meta.wcs.get_transform('slit_frame', 'msa_frame')
            is_fs = [s.quadrant == 5 for s in slit2msa.slits]
            if not any(is_fs):
                log.info("Masking the fixed slit region for MOS data.")
                mask[..., NRS_FS_REGION[0]:NRS_FS_REGION[1], :] = False
            else:
                log.info("Fixed slits found in MSA definition; "
                         "not masking the fixed slit region for MOS data.")

    # Mask outliers and signal using sigma clipping stats.
    # For 3D data, loop over each integration separately.
    if image_model.data.ndim == 3:
        for i in range(image_model.data.shape[0]):
            clip_to_background(
                image_model.data[i], mask[i],
                sigma_upper=n_sigma, fit_histogram=fit_histogram, verbose=True)
    else:
        clip_to_background(
            image_model.data, mask,
            sigma_upper=n_sigma, fit_histogram=fit_histogram, verbose=True)

    # Close the image model if needed
    if image_model is not input_model:
        image_model.close()
        del image_model

    # Reduce the mask to a single plane if needed
    if single_mask and mask.ndim == 3:
        mask = np.all(mask, axis=0)

    # Return the mask
    return mask


def background_level(image, mask, background_method='median'):
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
        background data is modeled with a median filter with a 5x5
        pixel kernel.  If None, the background value is 0.0.

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
            try:
                bkg = Background2D(
                    image, (34, 34), filter_size=(5, 5), mask=~mask,
                    sigma_clip=sigma_clip_for_bkg,
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


def median_clean(image, mask, slowaxis):
    """
    Fit and remove background noise via median values along the detector slow axis.

    Parameters
    ----------
    image : array-like of float
        The image to be cleaned.

    mask : array-like of bool
        The mask that indicates which pixels are to be used in fitting.
        True indicates a background pixel.

    slowaxis : int
        The detector slow readout direction.  Values expected are 1
        or 2, following the JWST datamodel definition (meta.subarray.slowaxis).

    Returns
    -------
    cleaned_image : array-like of float
        The cleaned image.
    """
    # Masked median along slow axis
    masked_image = np.ma.array(image, mask=~mask)
    stripes = np.ma.median(masked_image, axis=(slowaxis - 1), keepdims=True)
    stripes = np.ma.filled(stripes, fill_value=0.0)

    # Remove median stripes
    corrected_image = image - stripes

    return corrected_image


def do_correction(input_model, fit_method, background_method,
                  mask_spectral_regions, n_sigma, fit_histogram,
                  single_mask, save_mask, user_mask):
    """
    Apply the 1/f noise correction.

    Parameters
    ----------
    input_model : `~jwst.datamodel.JwstDataModel`
        Science data to be corrected.

    fit_method : {'fft', 'median'}
        The algorithm to use to fit background noise.

    background_method : {'median', 'model', None}
        If 'median', the preliminary background to remove and restore
        is a simple median of the background data.  If 'model', the
        background data is modeled with a median filter with a 5x5
        pixel kernel.  If None, the background value is 0.0.

    mask_spectral_regions : bool
        Mask slit/slice regions defined in WCS. Implemented only for
        NIRSpec science data modes.

    n_sigma : float
        N-sigma rejection level for finding outliers.

    fit_histogram : bool, optional
        If set, the 'sigma' used with `n_sigma` for clipping outliers
        is derived from a Gaussian fit to a histogram of values.
        Otherwise, a simple iterative sigma clipping is performed.

    single_mask : bool
        If set, a single mask will be created, regardless of
        the number of input integrations. Otherwise, the mask will
        be a 3D cube, with one plane for each integration.

    save_mask : bool
        Switch to indicate whether the mask should be saved.

    user_mask : str or None
        Path to user-supplied mask image.

    Returns
    -------
    output_model : `~jwst.datamodel.JwstDataModel`
        Corrected data.

    mask_model : `~jwst.datamodel.JwstDataModel`
        Pixel mask to be saved or None.
    """
    # Track the completion status, for various failure conditions
    status = 'SKIPPED'

    detector = input_model.meta.instrument.detector.upper()
    slowaxis = np.abs(input_model.meta.subarray.slowaxis)
    subarray = input_model.meta.subarray.name.upper()
    exp_type = input_model.meta.exposure.type
    log.info(f'Input exposure type is {exp_type}, detector={detector}')

    # Check for a valid input that we can work on
    nsclean_allowed = ['NRS_MSASPEC', 'NRS_IFU', 'NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']
    if fit_method == 'fft':
        message = None
        if exp_type not in nsclean_allowed:
            message = f"Fit method 'fft' cannot be applied to exp_type {exp_type}"
        elif subarray == 'ALLSLITS':
            message = f"Fit method 'fft' cannot be applied to subarray {subarray}"
        if message is not None:
            log.warning(message)
            log.warning("Step will be skipped")
            return input_model, None, status

    output_model = input_model.copy()

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
            input_model,
            mask_spectral_regions=mask_spectral_regions,
            n_sigma=n_sigma, fit_histogram=fit_histogram,
            single_mask=single_mask)

        # Store the mask image in a model, if requested
        if save_mask:
            if len(background_mask.shape) == 3:
                mask_model = datamodels.CubeModel(data=background_mask)
            else:
                mask_model = datamodels.ImageModel(data=background_mask)
        else:
            mask_model = None

    log.info(f"Cleaning image {input_model.meta.filename}")

    # Setup for handling 2D, 3D, or 4D inputs
    image_shape = input_model.data.shape[-2:]
    ndim = input_model.data.ndim
    if ndim == 2:
        nints = 1
        ngroups = 1
    else:
        nints = input_model.data.shape[0]
        if ndim == 3:
            ngroups = 1
        else:
            ngroups = input_model.data.shape[1] - 1

        # Check for 3D mask
        if background_mask.ndim == 2:
            if nints > 1:
                log.info("Data has multiple integrations, but mask is 2D: "
                         "the same mask will be used for all integrations.")
        elif background_mask.shape[0] != nints:
            log.warning("Mask does not match data shape. Step will be skipped.")
            return output_model, None, status

    # Check that mask matches 2D data shape
    if background_mask.shape[-2:] != image_shape:
        log.warning("Mask does not match data shape. Step will be skipped.")
        return output_model, None, status

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

            # Fit and remove a background level
            if str(background_method).lower() == 'none':
                background = 0.0
                bkg_sub = image
            else:
                background = background_level(
                    image, mask, background_method=background_method)
                log.debug(f'Background level: {np.nanmedian(background):.5g}')
                bkg_sub = image - background

                # Flag more signal in the background subtracted image,
                # with sigma set by the lower half of the distribution only
                clip_to_background(
                    bkg_sub, mask, sigma_lower=n_sigma,
                    sigma_upper=n_sigma, lower_half_only=True)

                # TODO: add segmentation masking?

            if fit_method == 'fft':
                if bkg_sub.shape == (2048, 2048):
                    # Clean a full-frame image
                    cleaned_image = fft_clean_full_frame(bkg_sub, mask, detector)
                else:
                    # Clean a subarray image
                    cleaned_image = fft_clean_subarray(bkg_sub, mask, detector)
            else:
                cleaned_image = median_clean(bkg_sub, mask, slowaxis)

            # Check for failure
            if cleaned_image is None:
                log.error(f'Cleaning failed for integration {i + 1}, group {j + 1}')
                # re-copy input to make sure any partial changes are thrown away
                output_model.close()
                del output_model
                output_model = input_model.copy()
                return output_model, None, status
            else:
                # Restore the background level
                cleaned_image += background

                # Restore NaNs in cleaned image
                cleaned_image[nan_pix] = np.nan

                # Store the cleaned image in the output model
                if ndim == 2:
                    output_model.data = cleaned_image
                elif ndim == 3:
                    output_model.data[i] = cleaned_image
                else:
                    # add the cleaned data diff to the previously cleaned group,
                    # rather than the noisy input group
                    output_model.data[i, j+1] = output_model.data[i, j] + cleaned_image

    # Set completion status
    status = 'COMPLETE'

    return output_model, mask_model, status
