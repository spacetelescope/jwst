import logging

import gwcs
import numpy as np
from gwcs.utils import to_index
from stdatamodels.jwst.datamodels import dqflags

from jwst import datamodels
from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst.clean_flicker_noise.background_level import background_level, clip_to_background
from jwst.clean_flicker_noise.nsclean import NSClean, NSCleanSubarray
from jwst.clean_flicker_noise.tso_median_image import make_median_image
from jwst.extract_1d.soss_extract import pastasoss
from jwst.flatfield import FlatFieldStep
from jwst.lib.basic_utils import disable_logging
from jwst.lib.reffile_utils import get_subarray_model, ref_matches_sci
from jwst.msaflagopen import MSAFlagOpenStep
from jwst.ramp_fitting import RampFitStep

log = logging.getLogger(__name__)

# Fixed slit region to mask, for NIRSpec MOS and IFU data
# Values are y start and stop indices, for the edges of the
# region to mask.
NRS_FS_REGION = [922, 1116]

__all__ = [
    "read_flat_file",
    "make_rate",
    "post_process_rate",
    "mask_ifu_slices",
    "mask_slits",
    "create_mask",
    "fft_clean_full_frame",
    "fft_clean_subarray",
    "median_clean",
    "make_intermediate_model",
    "do_correction",
]


def read_flat_file(input_model, flat_filename):
    """
    Read flat data from an input file path.

    Flat data is assumed to be full frame.  Subarrays matching the input
    data are extracted as needed.

    Only the flat image is returned: error and DQ arrays are ignored.
    Any zeros or NaNs in the flat image are set to a smoothed local average
    value (via ``background_level``, with ``background_method = 'model'``) before
    returning, to avoid impacting the background and noise fits near
    missing flat data.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input data.
    flat_filename : str
        File path for a full-frame flat image.

    Returns
    -------
    flat_data : ndarray of float
        A 2D flat image array matching the input data.
    """
    if flat_filename is None:
        return None

    # Open the provided flat as FlatModel
    log.debug("Dividing by flat data prior to fitting")
    flat = datamodels.FlatModel(flat_filename)

    # Extract subarray from reference data, if necessary
    if ref_matches_sci(input_model, flat):
        flat_data = flat.data
    else:
        log.debug("Extracting matching subarray from flat")
        sub_flat = get_subarray_model(input_model, flat)
        flat_data = sub_flat.data
        sub_flat.close()
    flat.close()

    # Set any zeros or non-finite values in the flat data to a smoothed local value
    bad_data = (flat_data == 0) | ~np.isfinite(flat_data)
    if np.any(bad_data):
        smoothed_flat = background_level(flat_data, ~bad_data, background_method="model")
        if np.isscalar(smoothed_flat):
            # 2D model failed, median value returned instead
            flat_data[bad_data] = smoothed_flat
        else:
            flat_data[bad_data] = smoothed_flat[bad_data]

    return flat_data


def make_rate(input_model, input_dir="", return_cube=False):
    """
    Make a rate model from a ramp model.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel`
        Input ramp model.
    input_dir : str
        Path to the input directory.  Used by sub-steps (e.g., ``assign_wcs``
        for NIRSpec MOS data) to find auxiliary data.
    return_cube : bool, optional
        If set, a `~stdatamodels.jwst.datamodels.CubeModel` will be returned, with a separate
        rate for each integration.  Otherwise, an
        `~stdatamodels.jwst.datamodels.ImageModel` is returned
        with the combined rate for the full observation.

    Returns
    -------
    rate_model : `~stdatamodels.jwst.datamodels.ImageModel` or \
                 `~stdatamodels.jwst.datamodels.CubeModel`
        The rate or rateints model.
    """
    # Call the ramp fit step on a copy of the input
    # Use software default values for parameters

    log.info("Creating draft rate file for scene masking")
    with disable_logging(level=logging.WARNING):
        step = RampFitStep()
        step.input_dir = input_dir

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


def post_process_rate(
    input_model, input_dir="", assign_wcs=False, msaflagopen=False, flat_dq=False
):
    """
    Perform additional processing for the input rate model, as needed.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.ImageModel` or \
                  `~stdatamodels.jwst.datamodels.CubeModel`
        Input rate model.
    input_dir : str
        Path to the input directory.  Used by sub-steps (e.g., ``assign_wcs``
        for NIRSpec MOS data) to find auxiliary data.
    assign_wcs : bool, optional
        If set and the input does not already have a WCS assigned,
        the ``assign_wcs`` step will be called on the rate model.
    msaflagopen : bool, optional
        If set, the ``msaflagopen`` step will be called on the rate model.
        If a WCS is not already present, ``assign_wcs`` will be called first.
    flat_dq : bool, optional
        If set, the ``flat_field`` step will be run on the input model. DQ
        flags are retrieved from the output and added to the input
        model's DQ array. The rate data is not modified.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.ImageModel` or \
                   `~stdatamodels.jwst.datamodels.CubeModel`
        The updated model.
    """
    output_model = input_model

    # If needed, assign a WCS
    if (assign_wcs or msaflagopen) and getattr(output_model.meta, "wcs", None) is None:
        log.info("Assigning a WCS for scene masking")
        with disable_logging(level=logging.WARNING):
            step = AssignWcsStep()
            step.input_dir = input_dir
            output_model = step.run(output_model)

    # If needed, flag open MSA shutters
    if msaflagopen:
        log.info("Flagging failed-open MSA shutters for scene masking")
        with disable_logging(level=logging.WARNING):
            step = MSAFlagOpenStep()
            step.input_dir = input_dir
            output_model = step.run(output_model)

    # If needed, draft a flat correction to retrieve non-science areas
    if flat_dq:
        log.info("Retrieving flat DQ values for scene masking")
        with disable_logging(level=logging.WARNING):
            step = FlatFieldStep()
            step.input_dir = input_dir
            flat_corrected_model = step.run(output_model)

        # Copy out the flat DQ plane, leave the data as is
        output_model.dq |= flat_corrected_model.dq
        flat_corrected_model.close()
        del flat_corrected_model

    return output_model


def mask_ifu_slices(input_model, mask):
    """
    Flag pixels within NIRSpec IFU slices.

    Find pixels located within IFU slices, according to the WCS,
    and flag them in the mask, so that they do not get used.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Science data model
    mask : ndarray of bool
        2D input mask that will be updated. `True` indicates background
        pixels to be used. IFU slice regions will be set to `False`.

    Returns
    -------
    mask : ndarray of bool
        2D output mask with additional flags for science pixels
    """
    log.info("Finding slice pixels for an IFU image")

    # Initialize global DQ map to all zero (OK to use)
    dqmap = np.zeros_like(input_model.dq)

    # Get the wcs objects for all IFU slices
    list_of_wcs = nirspec.nrs_ifu_wcs(input_model)

    # Loop over the IFU slices, finding the valid region for each
    for ifu_wcs in list_of_wcs:
        # Construct array indexes for pixels in this slice
        x, y = gwcs.wcstools.grid_from_bounding_box(ifu_wcs.bounding_box, step=(1, 1), center=True)

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
    Flag pixels within science regions for NIRSpec slit modes.

    Find pixels located within MOS or fixed slit footprints
    and flag them in the mask, so that they do not get used.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Science data model.
    mask : ndarray of bool
        2D input mask that will be updated. `True` indicates background
        pixels to be used. Slit regions will be set to `False`.

    Returns
    -------
    mask : ndarray of bool
        2D output mask with additional flags for slit pixels
    """
    log.info("Finding slit/slitlet pixels")

    # Get the slits from the WCS object
    slits = input_model.meta.wcs.get_transform("gwa", "slit_frame").slit_ids

    # Loop over the slits, marking all the pixels within each bounding
    # box as False (do not use) in the mask.
    # Note that for 3D masks (TSO mode), all planes will be set to the same value.
    for slit in slits:
        bbox = input_model.meta.wcs.bounding_box[slit]
        xlo, xhi = to_index(bbox[0])
        ylo, yhi = to_index(bbox[1])
        mask[..., ylo:yhi, xlo:xhi] = False

    return mask


def mask_soss_traces(input_model, mask, soss_refmodel=None, halfwidth=16):
    """
    Flag pixels within science regions for NIRISS SOSS.

    Get trace positions (dependent on pupil wheel position)
    and mask them with basic uniform width wings.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Science data model.
    mask : ndarray of bool
        2D input mask that will be updated. True indicates background
        pixels to be used. Trace regions will be set to False.
    soss_refmodel : `~stdatamodels.jwst.datamodels.PastasossModel`, optional
        PASTASOSS reference file model. If None, a default model
        will be retrieved.
    halfwidth : int, optional
        Width above and below the center of the trace to mask.

    Returns
    -------
    mask : ndarray of bool
        2D output mask with additional flags for trace pixels.
    """
    pwcpos = input_model.meta.instrument.pupil_position
    subarray = input_model.meta.subarray.name
    wavemaps, traces = pastasoss.get_soss_wavemaps(
        pwcpos, subarray=subarray, refmodel=soss_refmodel, spectraces=True
    )
    ny, nx = input_model.data.shape[-2:]
    yidx, xidx = np.mgrid[:ny, :nx]
    for xpts, ypts, _ in traces:
        # Sort by x
        order = np.argsort(xpts)
        xpts = xpts[order]
        ypts = ypts[order]

        # Interpolate trace points onto full array width
        x_full = np.arange(nx)
        y_full = np.interp(x_full, xpts, ypts)

        # Mask pixels within halfwidth above and below the trace
        y_low = np.floor(y_full - halfwidth)
        y_high = np.ceil(y_full + halfwidth)
        mask[..., (yidx >= y_low) & (yidx <= y_high)] = False

    return mask


def create_mask(
    input_model,
    mask_science_regions=False,
    n_sigma=2.0,
    fit_histogram=False,
    single_mask=False,
    soss_refmodel=None,
):
    """
    Create a mask identifying background pixels.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Science data model, containing rate data with all necessary
        pre-processing already performed.
    mask_science_regions : bool, optional
        For NIRSpec, mask regions of the image defined by WCS bounding
        boxes for slits/slices, as well as any regions known to be
        affected by failed-open MSA shutters. This requires the
        ``assign_wcs`` and ``msaflagopen`` steps to have been run on the
        input_model.
        For MIRI imaging, mask regions of the detector not used for science.
        This requires that DO_NOT_USE flags are set in the DQ array
        for the input_model.
        For NIRISS SOSS, use the ``soss_refmodel`` to mask the trace locations.
    n_sigma : float, optional
        Sigma threshold for masking outliers. Set to 0 to skip outlier rejection.
    fit_histogram : bool, optional
        If set, the 'sigma' used with ``n_sigma`` for clipping outliers
        is derived from a Gaussian fit to a histogram of values.
        Otherwise, a simple iterative sigma clipping is performed.
    single_mask : bool, optional
        If set, a single mask will be created, regardless of
        the number of input integrations. Otherwise, the mask will
        be a 3D cube, with one plane for each integration.
    soss_refmodel : `~stdatamodels.jwst.datamodels.PastasossModel`, optional
        If ``mask_science_regions`` is True and the input exposure
        is NIS_SOSS, will be used to mask the SOSS traces.

    Returns
    -------
    mask : ndarray of bool
        2D or 3D image mask
    """
    exptype = input_model.meta.exposure.type.lower()

    # Initialize mask to all True. Subsequent operations will mask
    # out pixels that contain signal.
    mask = np.full(input_model.dq.shape, True)

    # If IFU, mask all pixels contained in the IFU slices
    if exptype == "nrs_ifu" and mask_science_regions:
        mask = mask_ifu_slices(input_model, mask)

    # If MOS or FS, mask all pixels affected by open slitlets
    if exptype in ["nrs_fixedslit", "nrs_brightobj", "nrs_msaspec"] and mask_science_regions:
        mask = mask_slits(input_model, mask)

    # If NIRSpec, mask pixels affected by failed-open shutters
    if mask_science_regions and exptype.startswith("nrs"):
        open_pix = input_model.dq & dqflags.pixel["MSA_FAILED_OPEN"]
        mask[open_pix > 0] = False

    # If MIRI imaging, mask the non-science regions:
    # they contain irrelevant emission
    if mask_science_regions and exptype in ["mir_image"]:
        non_science = (
            input_model.dq & (dqflags.pixel["DO_NOT_USE"] | dqflags.pixel["NON_SCIENCE"])
        ) > 0
        mask[non_science] = False

    # If NIRISS SOSS, mask the traces from the reference model
    if mask_science_regions and exptype == "nis_soss":
        mask = mask_soss_traces(input_model, mask, soss_refmodel=soss_refmodel)

    # Mask any NaN pixels or exactly zero value pixels
    no_data_pix = np.isnan(input_model.data) | (input_model.data == 0)
    mask[no_data_pix] = False

    # If IFU or MOS, mask the fixed-slit area of the image; uses hardwired indexes
    if mask_science_regions:
        if exptype == "nrs_ifu":
            log.info("Masking the fixed slit region for IFU data.")
            mask[..., NRS_FS_REGION[0] : NRS_FS_REGION[1], :] = False
        elif exptype == "nrs_msaspec":
            # check for any slits defined in the fixed slit quadrant:
            # if there is nothing there of interest, mask the whole FS region
            try:
                slits = input_model.meta.wcs.get_transform("gwa", "slit_frame").slits
                is_fs = [s.quadrant == 5 for s in slits]
            except (AttributeError, ValueError, TypeError):
                log.warning("Open slits not found in input WCS.")
                is_fs = [False]
            if not any(is_fs):
                log.info("Masking the fixed slit region for MOS data.")
                mask[..., NRS_FS_REGION[0] : NRS_FS_REGION[1], :] = False
            else:
                log.info(
                    "Fixed slits found in MSA definition; "
                    "not masking the fixed slit region for MOS data."
                )

    # Mask outliers and signal using sigma clipping stats.
    # For 3D data, loop over each integration separately.
    if n_sigma > 0:
        if input_model.data.ndim == 3:
            for i in range(input_model.data.shape[0]):
                clip_to_background(
                    input_model.data[i],
                    mask[i],
                    sigma_upper=n_sigma,
                    fit_histogram=fit_histogram,
                    verbose=True,
                )
        else:
            clip_to_background(
                input_model.data,
                mask,
                sigma_upper=n_sigma,
                fit_histogram=fit_histogram,
                verbose=True,
            )

    # Reduce the mask to a single plane if needed
    if single_mask and mask.ndim == 3:
        mask = np.all(mask, axis=0)

    # Return the mask
    return mask


def fft_clean_full_frame(image, mask, detector):
    """
    Fit and remove background noise in frequency space for a full-frame image.

    Parameters
    ----------
    image : ndarray of float
        The image to be cleaned.
    mask : ndarray of bool
        The mask that indicates which pixels are to be used in fitting.
    detector : str
        The name of the detector from which the data originate.

    Returns
    -------
    cleaned_image : ndarray of float
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


def fft_clean_subarray(
    image,
    mask,
    detector,
    npix_iter=512,
    fc=(1061, 1211, 49943, 49957),
    exclude_outliers=True,
    sigrej=4,
    minfrac=0.05,
):
    """
    Fit and remove background noise in frequency space for a subarray image.

    Parameters
    ----------
    image : ndarray of float
        The image to be cleaned.
    mask : ndarray of bool
        The mask that indicates which pixels are to be used in fitting.
    detector : str
        The name of the detector from which the data originate.
    npix_iter : int
        Number of pixels to process simultaneously.  Default 512.  Should
        be at least a few hundred to access sub-kHz frequencies in areas
        where most pixels are available for fitting.  Previous default
        behavior corresponds to ``npix_iter`` of infinity.
    fc : tuple
        Apodizing filter definition. These parameters are tunable. The
        defaults happen to work well for NIRSpec BOTS exposures:

        1. Unity gain for ``f < fc[0]``
        2. Cosine roll-off from ``fc[0]`` to ``fc[1]``
        3. Zero gain from ``fc[1]`` to ``fc[2]``
        4. Cosine roll-on from ``fc[2]`` to ``fc[3]``

        Default ``(1061, 1211, 49943, 49957)``
    exclude_outliers : bool
        Find and mask outliers in the fit? Default: `True`
    sigrej : float
        Number of sigma to clip when identifying outliers. Default: 4.
    minfrac : float
        Minimum fraction of pixels locally available in the mask in
        order to attempt a correction. Default: 0.05 (i.e., 5%)

    Returns
    -------
    cleaned_image : ndarray of float
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

    # We must do the masking of discrepant pixels here: it just
    # doesn't work if we wait and do it in the cleaner.  This is
    # basically copied from lib.py.  Use a robust estimator for
    # standard deviation, then exclude discrepant pixels and their
    # four nearest neighbors from the fit.

    if exclude_outliers:
        med = np.median(image[mask])
        std = 1.4825 * np.median(np.abs((image - med)[mask]))
        outlier = mask & (np.abs(image - med) > sigrej * std)

        mask = mask & (~outlier)

        # also get four nearest neighbors of flagged pixels
        mask[1:] = mask[1:] & (~outlier[:-1])
        mask[:-1] = mask[:-1] & (~outlier[1:])
        mask[:, 1:] = mask[:, 1:] & (~outlier[:, :-1])
        mask[:, :-1] = mask[:, :-1] & (~outlier[:, 1:])

    # Used to determine the fitting intervals along the slow scan
    # direction.  Prepend a zero so that sum_mask[i] is equal
    # to np.sum(mask[:i], axis=1).

    sum_mask = np.array([0] + list(np.cumsum(np.sum(mask, axis=1))))

    # i1 will be the first row with a nonzero element in the mask
    # imax will be the last row with a nonzero element in the mask

    nonzero_mask_element = np.sum(mask, axis=1) > 0

    if np.sum(nonzero_mask_element) == 0:
        log.warning("No good pixels in mask; step will be skipped")
        return None

    i1 = np.amin(np.arange(mask.shape[0])[nonzero_mask_element])
    imax = np.amax(np.arange(mask.shape[0])[nonzero_mask_element])

    i1_vals = []
    di_list = []
    models = []
    while i1 <= imax:
        # Want npix_iter available pixels in this section.  If
        # there are fewer than 1.5*npix_iter available pixels in
        # the rest of the image, just go to the end.
        k = 0
        for k in range(i1 + 1, imax + 2):
            if (
                sum_mask[k] - sum_mask[i1] > npix_iter
                and sum_mask[-1] - sum_mask[i1] > 1.5 * npix_iter
            ):
                break

        di = k - i1

        i1_vals += [i1]
        di_list += [di]

        # Fit this section only if at least minpct% of the pixels
        # are available for finding the background.  Don't flag
        # outliers section-by-section; we have to do that earlier
        # over the full array to get reliable values for the mean
        # and standard deviation.

        if np.mean(mask[i1 : i1 + di]) > minfrac:
            cleaner = NSCleanSubarray(
                image[i1 : i1 + di], mask[i1 : i1 + di], fc=fc, exclude_outliers=False
            )
            try:
                models += [cleaner.clean(return_model=True)]
            except np.linalg.LinAlgError:
                log.warning("Error cleaning image; step will be skipped")
                return None
        else:
            log.warning(
                "Insufficient reference pixels for NSClean around "
                f"row {i1}; no correction will be made here."
            )
            models += [np.zeros(image[i1 : i1 + di].shape)]

        # If we have reached the end of the array, we are finished.
        if k == imax + 1:
            break

        # Step forward by half an interval so that we have
        # overlapping fitting regions.

        i1 += max(int(np.round(di / 2)), 1)

    model = np.zeros(image.shape)
    tot_wgt = np.zeros(image.shape)

    # When we combine different corrections computed over
    # different intervals, each one the highest weight towards the
    # center of its interval and less weight towards the edge.
    # Use nonzero weights everywhere so that if only one
    # correction is available it gets unit weight when we
    # normalize.

    for i in range(len(models)):
        wgt = 1.001 - np.abs(np.linspace(-1, 1, di_list[i]))[:, np.newaxis]
        model[i1_vals[i] : i1_vals[i] + di_list[i]] += wgt * models[i]
        tot_wgt[i1_vals[i] : i1_vals[i] + di_list[i]] += wgt

    # don't divide by zero
    tot_wgt[model == 0] = 1
    model /= tot_wgt
    cleaned_image = image - model

    # Restore the cleaned image to the science frame
    if detector == "NRS1":
        cleaned_image = cleaned_image.transpose()
    else:
        cleaned_image = cleaned_image[::-1].transpose()

    return cleaned_image


def median_clean(image, mask, axis_to_correct, fit_by_channel=False):
    """
    Fit and remove background noise via median values along one image axis.

    Parameters
    ----------
    image : ndarray of float
        The image to be cleaned.
    mask : ndarray of bool
        The mask that indicates which pixels are to be used in fitting.
        True indicates a background pixel.
    axis_to_correct : int
        For NIR detectors, the axis to correct should be the detector slow
        readout direction.  Values expected are 1 or 2, following the JWST
        datamodel definition (``meta.subarray.slowaxis``).  For MIRI, flicker
        noise appears along the vertical direction, so ``axis_to_correct``
        should be set to 1 (median along the y-axis).
    fit_by_channel : bool, optional
        If set, flicker noise is fit independently for each detector channel.
        Ignored for MIRI and for subarray data.

    Returns
    -------
    cleaned_image : ndarray of float
        The cleaned image.
    """
    # Masked median along slow axis
    masked_image = np.ma.array(image, mask=~mask)
    corrected_image = image.copy()
    array_axis = axis_to_correct - 1

    # If desired, take the median over each channel separately.

    # Full frame for the NIR detectors is 2048 x 2048, and there
    # are 4 channels per detector, divided along the slow axis.
    # The fit_by_channel option should only be set for full frame
    # NIR data, but make sure the slow axis has the right size
    # here anyway.
    if fit_by_channel and image.shape[array_axis] == 2048:
        n_output = 4
        channel_size = 512
    else:
        # For MIRI data and subarrays, the whole image should
        # be considered at once, as a single "channel"
        n_output = 1
        channel_size = image.shape[array_axis]

    # Compute stripes from the median over the background data
    # in the channel
    cstart = 0
    cstop = channel_size
    for _channel in range(n_output):
        if array_axis == 1:
            channel_image = masked_image[:, cstart:cstop]
        else:
            channel_image = masked_image[cstart:cstop, :]
        stripes = np.ma.median(channel_image, axis=array_axis, keepdims=True)
        stripes = np.ma.filled(stripes, fill_value=0.0)

        # Remove median stripes
        if array_axis == 1:
            corrected_image[:, cstart:cstop] = image[:, cstart:cstop] - stripes
        else:
            corrected_image[cstart:cstop, :] = image[cstart:cstop, :] - stripes

        cstart += channel_size
        cstop += channel_size

    return corrected_image


def make_intermediate_model(input_model, intermediate_data):
    """
    Make a data model to contain intermediate outputs.

    The output model type depends on the shape of the input
    intermediate data.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input data.
    intermediate_data : ndarray
        The intermediate data to save.

    Returns
    -------
    intermediate_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        A model containing only the intermediate data and top-level
        metadata matching the input.
    """
    if intermediate_data.ndim == 4:
        intermediate_model = datamodels.RampModel(data=intermediate_data)
    elif intermediate_data.ndim == 3:
        intermediate_model = datamodels.CubeModel(data=intermediate_data)
    else:
        intermediate_model = datamodels.ImageModel(data=intermediate_data)

    # Copy metadata from input model
    intermediate_model.update(input_model)
    return intermediate_model


def _check_input(exp_type, fit_method):
    """
    Check for valid input data and options.

    Parameters
    ----------
    exp_type : str
        Exposure type for the input.
    fit_method : str
        Noise fitting method.

    Returns
    -------
    valid : bool
        True if the input is valid.
    """
    message = None
    miri_allowed = ["MIR_IMAGE"]
    if exp_type.startswith("MIR") and exp_type not in miri_allowed:
        message = f"EXP_TYPE {exp_type} is not supported"

    nsclean_allowed = ["NRS_MSASPEC", "NRS_IFU", "NRS_FIXEDSLIT", "NRS_BRIGHTOBJ"]
    if fit_method == "fft":
        if exp_type not in nsclean_allowed:
            message = f"Fit method 'fft' cannot be applied to exp_type {exp_type}"

    if message is not None:
        log.warning(message)
        log.warning("Step will be skipped")
        return False

    return True


def _standardize_parameters(
    exp_type, subarray, slowaxis, background_method, fit_by_channel, single_mask
):
    """
    Standardize input parameters.

    Check input parameters against input exposure type and assemble
    values needed for subsequent correction.

    Parameters
    ----------
    exp_type : str
        Exposure type for the input data.
    subarray : str
        Subarray name for the input data.
    slowaxis : int
        Detector slow axis.
    background_method : str
        Input option for background method.
    fit_by_channel : bool
        Input option to fit noise by channel.
    single_mask : bool
        Input option to make one mask for all integrations.

    Returns
    -------
    axis_to_correct : int
        Axis along which flicker noise appears for the input
        exposure type.
    background_method : str or None
        Standardized parameter value for background correction.
    fit_by_channel : bool
        Standardized parameter value for the input subarray and
        exposure type.
    single_mask : bool
        Standardized parameter value for the background method.
    fc : tuple of int
        Frequency cutoff values for use with FFT correction,
        by input subarray.
    """
    # Get axis to correct, by instrument
    if exp_type.startswith("MIR"):
        # MIRI doesn't have 1/f-noise, but it does have a vertical flickering.
        # Set the axis for median correction to the y-axis.
        axis_to_correct = 1
    else:
        # For NIR detectors, the axis for median correction is the
        # detector slowaxis.
        axis_to_correct = abs(slowaxis)

    # Standardize background arguments
    if str(background_method).lower() == "none":
        background_method = None
    if str(background_method).lower() == "median_image" and single_mask:
        log.warning("The median_image background method requires draft rateints files.")
        log.warning("Setting single_mask to False.")
        single_mask = False

    # Check for fit_by_channel argument, and use only if data is full frame
    if fit_by_channel and (subarray != "FULL" or exp_type.startswith("MIR")):
        log.warning("Fit by channel can only be used for full-frame NIR data.")
        log.warning("Setting fit_by_channel to False.")
        fit_by_channel = False

    # For a fft correction, ALLSLITS exposures need different
    # ranges of 1/f frequencies.  Be less aggressive with
    # fitting higher frequencies in ALLSLITS mode.
    if subarray == "ALLSLITS":
        fc = (150, 200, 49943, 49957)
    else:
        fc = (1061, 1211, 49943, 49957)

    return axis_to_correct, background_method, fit_by_channel, single_mask, fc


def _make_processed_rate_image(
    input_model, single_mask, input_dir, exp_type, mask_science_regions, flat, background_method
):
    """
    Make a draft rate image and postprocess if needed.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input data.
    single_mask : bool
        If set, a single scene mask is desired, so create
        a single rate image.  Otherwise, create a rateints
        cube, with independent rate information for each
        integration.
    input_dir : str
        Path to the input directory.  Used by sub-steps to find auxiliary
        data.
    exp_type : str
        Exposure type for the input data, used to determine
        which kinds of postprocessing is necessary.
    mask_science_regions : bool
        If set, science regions should be masked, so run ``assign_wcs``
        and ``msaflagopen`` for NIRSpec and ``flatfield`` for MIRI
        after creating the draft rate image.
    flat : ndarray of float or None
        If not None, the draft rate will be divided by the flat array
        before returning. The provided flat must match the rate shape.
    background_method : str
        If 'median_image', ``assign_wcs`` needs to be run for most
        spectrosopic modes.

    Returns
    -------
    image_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The processed rate image or cube.
    """
    if isinstance(input_model, datamodels.RampModel):
        image_model = make_rate(input_model, return_cube=(not single_mask), input_dir=input_dir)
    else:
        # input is already a rate file
        image_model = input_model

    # If needed, assign a WCS to the rate file,
    # flag open MSA shutters, or retrieve flat DQ flags
    assign_wcs = exp_type.startswith("NRS") or (
        background_method == "median_image" and exp_type != "NIS_SOSS"
    )
    flag_open = exp_type in ["NRS_IFU", "NRS_MSASPEC"]
    flat_dq = exp_type in ["MIR_IMAGE"]
    if mask_science_regions or background_method == "median_image":
        image_model = post_process_rate(
            image_model,
            assign_wcs=assign_wcs,
            msaflagopen=flag_open,
            flat_dq=flat_dq,
            input_dir=input_dir,
        )

    # Divide by the flat if provided
    if flat is not None:
        # Make a copy if necessary so we don't change the input rate file
        if image_model is input_model:
            image_model = input_model.copy()
        image_model.data /= flat

    return image_model


def _make_scene_mask(
    user_mask,
    image_model,
    mask_science_regions,
    n_sigma,
    fit_histogram,
    single_mask,
    save_mask,
    soss_refmodel=None,
):
    """
    Make a scene mask from user input or rate image.

    If provided, the user mask is opened as a datamodel and directly
    returned. Otherwise, the mask is generated from the rate data in
    ``image_model``.

    Parameters
    ----------
    user_mask : str or None
        Path to user-supplied mask image.
    image_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        A rate image or cube, processed as needed.
    mask_science_regions : bool
        For NIRSpec, mask regions of the image defined by WCS bounding
        boxes for slits/slices, as well as any regions known to be
        affected by failed-open MSA shutters.  For MIRI imaging, mask
        regions of the detector not used for science.
    n_sigma : float
        N-sigma rejection level for finding outliers. Set to 0 to skip
        outlier rejection.
    fit_histogram : bool
        If set, the 'sigma' used with ``n_sigma`` for clipping outliers
        is derived from a Gaussian fit to a histogram of values.
        Otherwise, a simple iterative sigma clipping is performed.
    single_mask : bool
        If set, a single mask will be created, regardless of
        the number of input integrations. Otherwise, the mask will
        be a 3D cube, with one plane for each integration.
    save_mask : bool
        If set, a mask model is created and returned along with the mask
        array. If not, the ``mask_model`` returned is None.
    soss_refmodel : `~stdatamodels.jwst.datamodels.PastasossModel`, optional
        If ``mask_science_regions`` is True and the input exposure
        is NIS_SOSS, will be used to mask the SOSS traces.

    Returns
    -------
    background_mask : ndarray of bool
        Mask array, with True indicating background pixels, False
        indicating source pixels.
    mask_model : `~stdatamodels.jwst.datamodels.JwstDataModel` or None
        A datamodel containing the background mask, if ``save_mask``
        is `True`.
    """
    # Check for a user-supplied mask image. If provided, use it.
    if user_mask is not None:
        mask_model = datamodels.open(user_mask)
        background_mask = (mask_model.data.copy()).astype(np.bool_)
        mask_model.close()
        del mask_model
    else:
        # Create the pixel mask that will be used to indicate which pixels
        # to include in the 1/f noise measurements. Basically, we're setting
        # all illuminated pixels to False, so that they do not get used, and
        # setting all unilluminated pixels to True (so they DO get used).
        log.info("Creating mask")
        background_mask = create_mask(
            image_model,
            mask_science_regions=mask_science_regions,
            n_sigma=n_sigma,
            fit_histogram=fit_histogram,
            single_mask=single_mask,
            soss_refmodel=soss_refmodel,
        )

    # Store the mask image in a model, if requested
    if save_mask:
        mask_model = make_intermediate_model(image_model, background_mask)
    else:
        mask_model = None

    return background_mask, mask_model


def _check_data_shapes(input_model, background_mask):
    """
    Check data shape for input model and background mask.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input data model.
    background_mask : ndarray of bool
        The background mask.

    Returns
    -------
    mismatch : bool
        If `True`, the background mask does not match the data
        and the step should be skipped.
    ndim : int
        Number of dimensions in the input data.
    nints : int
        Number of integrations in the input data.
    ngroups : int
        Number of groups in the input data.
    """
    mismatch = True
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
                log.info(
                    "Data has multiple integrations, but mask is 2D: "
                    "the same mask will be used for all integrations."
                )
        elif background_mask.shape[0] != nints:
            log.warning("Mask does not match data shape. Step will be skipped.")
            return mismatch, ndim, nints, ngroups

    # Check that mask matches 2D data shape
    if background_mask.shape[-2:] != image_shape:
        log.warning("Mask does not match data shape. Step will be skipped.")
        return mismatch, ndim, nints, ngroups

    return False, ndim, nints, ngroups


def _clean_one_image(
    image,
    mask,
    background_method,
    background_box_size,
    n_sigma,
    fit_method,
    detector,
    fc,
    axis_to_correct,
    fit_by_channel,
    flat,
):
    """
    Clean an image by fitting and removing background noise.

    Parameters
    ----------
    image : ndarray of float
        The input image.
    mask : ndarray of bool
        The scene mask.  All non-background signal should be
        marked as `True`.
    background_method : str
        The method for fitting the background.
    background_box_size : tuple of int
        Background box size, used when ``background_method``
        is 'model'.
    n_sigma : float
        N-sigma rejection level for outliers.
    fit_method : str
        The method for fitting the noise after background
        signal is removed.
    detector : str
        The detector name.
    fc : tuple of int
        Frequency cutoff values, used for ``fit_method = 'fft'``.
    axis_to_correct : int
        Axis along which noise appears, used for
        ``fit_method = 'median'``.
    fit_by_channel : bool
        Option to fit noise in detector channels separately,
        used for ``fit_method = 'median'``
    flat : ndarray of float or None
        If not None, the image is divided by the flat before fitting background
        and noise.  Flat data must match the shape of the image.

    Returns
    -------
    cleaned_image : ndarray of float
        The cleaned image.
    background : ndarray of float
        The background fit and removed prior to cleaning.
        Used for diagnostic purposes.
    success : bool
        `True` if cleaning proceeded as expected; `False` if
        cleaning failed and the step should be skipped.
    """
    success = True

    # Make sure data is float32
    image = image.astype(np.float32)

    # If provided, divide the image by the flat
    if flat is not None:
        image /= flat

    # Find and replace/mask NaNs
    nan_pix = np.isnan(image)
    image[nan_pix] = 0.0
    mask[nan_pix] = False

    # If there's no good data remaining in the mask,
    # skip this image
    if not np.any(mask):
        return None, None, success

    # Fit and remove a background level
    if str(background_method).lower() == "none":
        background = 0.0
        bkg_sub = image
    else:
        background = background_level(
            image,
            mask,
            background_method=background_method,
            background_box_size=background_box_size,
        )
        log.debug(f"Background level: {np.nanmedian(background):.5g}")
        bkg_sub = image - background

        # Flag more signal in the background subtracted image,
        # with sigma set by the lower half of the distribution only
        if n_sigma > 0:
            clip_to_background(
                bkg_sub, mask, sigma_lower=n_sigma, sigma_upper=n_sigma, lower_half_only=True
            )

    # Clean the noise
    if fit_method == "fft":
        if bkg_sub.shape == (2048, 2048):
            # Clean a full-frame image
            cleaned_image = fft_clean_full_frame(bkg_sub, mask, detector)
        else:
            # Clean a subarray image
            cleaned_image = fft_clean_subarray(bkg_sub, mask, detector, fc=fc)
    else:
        cleaned_image = median_clean(bkg_sub, mask, axis_to_correct, fit_by_channel=fit_by_channel)
    if cleaned_image is None:
        return None, None, False

    # Restore the background level
    cleaned_image += background

    # Restore the flat structure
    if flat is not None:
        cleaned_image *= flat

    # Restore NaNs in cleaned image
    cleaned_image[nan_pix] = np.nan

    return cleaned_image, background, success


def _mask_unusable(mask, dq):
    """
    Mask unusable data, according to DQ.

    Currently, JUMP and DO_NOT_USE flags are used.
    The mask is updated in place

    Parameters
    ----------
    mask : ndarray of bool
        Input mask, updated in place.
    dq : ndarray of int
        DQ flag array matching the mask shape.
    """
    dnu = (dq & dqflags.group["DO_NOT_USE"]) > 0
    mask[dnu] = False
    jump = (dq & dqflags.group["JUMP_DET"]) > 0
    mask[jump] = False


def do_correction(
    input_model,
    input_dir=None,
    fit_method="median",
    fit_by_channel=False,
    background_method="median",
    background_box_size=None,
    mask_science_regions=False,
    flat_filename=None,
    pastasoss_filename=None,
    n_sigma=2.0,
    fit_histogram=False,
    single_mask=True,
    user_mask=None,
    save_mask=False,
    save_background=False,
    save_noise=False,
):
    """
    Apply the 1/f noise correction.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Science data to be corrected. Updated in place.
    input_dir : str
        Path to the input directory.  Used by sub-steps (e.g., ``assign_wcs``
        for NIRSpec MOS data) to find auxiliary data.
    fit_method : {'fft', 'median'}, optional
        The algorithm to use to fit background noise.
    fit_by_channel : bool, optional
        If set, flicker noise is fit independently for each detector channel.
        Ignored for MIRI, for subarray data, and for ``fit_method = 'fft'``.
    background_method : {'median', 'model', 'median_image', None}, optional
        If 'median', the preliminary background to remove and restore
        is a simple median of the background data.  If 'model', the
        background data is fit with a low-resolution model via
        `~photutils.background.Background2D`.  If 'median_image' and the input
        has multiple integrations, a median image across integrations will
        be computed and subtracted as the background. If None, the background
        value is 0.0.
    background_box_size : tuple of int, optional
        Box size for the data grid used by `~photutils.background.Background2D` when
        ``background_method = 'model'``. For best results, use a box size
        that evenly divides the input image shape.
    mask_science_regions : bool, optional
        For NIRSpec, mask regions of the image defined by WCS bounding
        boxes for slits/slices, as well as any regions known to be
        affected by failed-open MSA shutters.  For MIRI imaging, mask
        regions of the detector not used for science.
    flat_filename : str, optional
        Path to a flat field image to apply to the data before fitting
        noise/background.
    pastasoss_filename : str, optional
        Path to a ``pastasoss`` reference file name. Used for NIS_SOSS only.
    n_sigma : float, optional
        N-sigma rejection level for finding outliers. Set to 0 to skip
        outlier rejection.
    fit_histogram : bool, optional
        If set, the 'sigma' used with ``n_sigma`` for clipping outliers
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
    output_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Corrected data.
    mask_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Pixel mask to be saved or None.
    background_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Background model to be saved or None.
    noise_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Background model to be saved or None.
    status : {'COMPLETE', 'SKIPPED'}
        Completion status.  If errors were encountered, status is set to 'SKIPPED'
        and the output data matches the input data.  Otherwise,
        status is set to 'COMPLETE'.
    """
    # Track the completion status, for various failure conditions
    status = "SKIPPED"

    detector = input_model.meta.instrument.detector.upper()
    subarray = input_model.meta.subarray.name.upper()
    exp_type = input_model.meta.exposure.type
    slowaxis = input_model.meta.subarray.slowaxis
    log.info(f"Input exposure type is {exp_type}, detector={detector}")

    # Check for a valid input that we can work on
    if not _check_input(exp_type, fit_method):
        return input_model, None, None, None, status

    # Get parameters needed for subsequent corrections, as appropriate
    # to the input data
    axis_to_correct, background_method, fit_by_channel, single_mask, fc = _standardize_parameters(
        exp_type, subarray, slowaxis, background_method, fit_by_channel, single_mask
    )

    # Read the flat and pastasoss files, if provided
    flat = read_flat_file(input_model, flat_filename)
    if pastasoss_filename is not None:
        soss_refmodel = datamodels.PastasossModel(pastasoss_filename)
    else:
        soss_refmodel = None

    # Make a rate file if needed
    if user_mask is None or background_method == "median_image":
        image_model = _make_processed_rate_image(
            input_model,
            single_mask,
            input_dir,
            exp_type,
            mask_science_regions,
            flat,
            background_method,
        )
    else:
        image_model = input_model

    # Make a mask model from the user input or the rate data
    background_mask, mask_model = _make_scene_mask(
        user_mask,
        image_model,
        mask_science_regions,
        n_sigma,
        fit_histogram,
        single_mask,
        save_mask,
        soss_refmodel=soss_refmodel,
    )

    # Check data shapes for 2D, 3D, or 4D inputs
    mismatch, ndim, nints, ngroups = _check_data_shapes(input_model, background_mask)
    if mismatch:
        return input_model, None, None, None, status

    # Make a background cube from a median image across integrations if needed
    if background_method == "median_image":
        if nints <= 1:
            log.warning(f"The median_image background method cannot be used with nints={nints}")
            log.warning("The step will be skipped.")
            return input_model, None, None, None, status
        try:
            median_image = make_median_image(input_model, image_model, soss_refmodel=soss_refmodel)
        except ValueError as err:
            log.warning("A median image could not be created.")
            log.warning(f"The error was: {err}")
            log.warning("The step will be skipped.")
            return input_model, None, None, None, status

        # The median image will be directly subtracted. Set background
        # method to None so any residual variable background levels are removed.
        background_method = None

        # Also add one to the ngroups for ramp data, since the first group
        # can now be processed
        if ndim > 3:
            ngroups += 1
    else:
        median_image = None

    # Close the draft rate model if created - it is no longer needed.
    if image_model is not input_model:
        image_model.close()
        del image_model

    # Make a background cube for saving, if desired
    if save_background:
        background_to_save = np.zeros_like(input_model.data)
    else:
        background_to_save = None

    # Keep a copy of the original input data
    input_data = input_model.data.copy()

    # Loop over integrations and groups (even if there's only 1)
    log.info(f"Cleaning image {input_model.meta.filename}")
    for i in range(nints):
        log.debug(f"Working on integration {i + 1} with ngroups {ngroups}")
        for j in range(ngroups):
            # Copy the scene mask, for further flagging
            if background_mask.ndim == 3:
                mask = background_mask[i].copy()
            else:
                mask = background_mask.copy()

            # Get the relevant image data
            if ndim == 2:
                image = input_data
            elif ndim == 3:
                if median_image is not None:
                    image = input_data[i] - median_image[i]
                else:
                    image = input_data[i]
            else:
                # Ramp data input
                if median_image is not None:
                    # subtract the median ramp
                    image = input_data[i, j] - median_image[i, j]
                else:
                    # subtract the current group from the next one
                    image = input_data[i, j + 1] - input_data[i, j]
                    dq = input_model.groupdq[i, j + 1]

                    # Mask any DNU and JUMP pixels
                    _mask_unusable(mask, dq)

            # Clean the image
            cleaned_image, background, success = _clean_one_image(
                image,
                mask,
                background_method,
                background_box_size,
                n_sigma,
                fit_method,
                detector,
                fc,
                axis_to_correct,
                fit_by_channel,
                flat,
            )

            if not success:
                # Cleaning failed for internal reasons - probably the
                # mask is not a good match to the data.
                log.error(f"Cleaning failed for integration {i + 1}, group {j + 1}")

                # Restore input data to make sure any partial changes
                # are thrown away
                input_model.data = input_data
                return input_model, None, None, None, status

            if cleaned_image is None:
                # Cleaning did not proceed because the image is bad:
                # leave it as is but continue correcting the rest
                log.warning(
                    f"No usable data in integration {i + 1}, group {j + 1}. "
                    f"Skipping correction for this image."
                )
                continue

            # Store the cleaned image in the output model
            if ndim == 2:
                input_model.data = cleaned_image
                if save_background:
                    background_to_save[:] = background
            elif ndim == 3:
                if median_image is not None:
                    # Add the median image back in
                    input_model.data[i] = median_image[i] + cleaned_image
                    if save_background:
                        background_to_save[i] = median_image[i]
                else:
                    input_model.data[i] = cleaned_image
                    if save_background:
                        background_to_save[i] = background
            else:
                if median_image is not None:
                    # Add the median ramp back in
                    input_model.data[i, j] = median_image[i, j] + cleaned_image
                    if save_background:
                        background_to_save[i, j] = median_image[i, j]
                else:
                    # Add the cleaned data diff to the previously cleaned group,
                    # rather than the noisy input group
                    input_model.data[i, j + 1] = input_model.data[i, j] + cleaned_image
                    if save_background:
                        background_to_save[i, j + 1] = background

    # Store the background image in a model, if requested
    if save_background:
        background_model = make_intermediate_model(input_model, background_to_save)
    else:
        background_model = None

    # Make a fit noise model for diagnostic purposes by
    # diffing the input and output models
    if save_noise:
        noise_data = input_model.data - input_data
        noise_model = make_intermediate_model(input_model, noise_data)
    else:
        noise_model = None

    # Set completion status
    status = "COMPLETE"

    return input_model, mask_model, background_model, noise_model, status
