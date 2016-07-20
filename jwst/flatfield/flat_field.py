#
#  Module for applying flat fielding
#

from __future__ import division

import numpy as np
import logging
from .. import datamodels
from .. datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def do_correction(input_model, flat_model, flat_suffix=None):
    """
    Short Summary
    -------------
    Flat-field a JWST data model using a flat-field model

    Parameters
    ----------
    input_model: JWST data model
        input science data model to be flat-fielded

    flat_model: JWST data model
        data model containing flat-field

    flat_suffix: str or None
        Filename suffix for optional output file to save flat field images.
        Note that this is only supported for MOS data.

    Returns
    -------
    output_model, interpolated_flats
        output_model is the data model for the flat-fielded science data
        interpolated_flats may be None, or it may be a MultiSlitModel
        containing the interpolated flat fields (MOS data only)
    """

    # Initialize the output model as a copy of the input
    output_model = input_model.copy()

    if isinstance(flat_model, datamodels.CubeFlatModel):
        interpolated_flats = do_MOS_flat_field(output_model,
                                               flat_model, flat_suffix)
    else:
        if flat_suffix is not None:
            log.warning("The flat_suffix parameter is not implemented"
                        " for this mode; will be ignored.")
        do_slit_flat_field(output_model, flat_model)
        interpolated_flats = None

    return (output_model, interpolated_flats)

#
# These functions are for slit (i.e. not MSA) flat fielding.
#
def do_slit_flat_field(output_model, flat_model):
    """
    Short Summary
    -------------
    Apply flat-fielding for the slit (not MSA) mode, updating the output model.

    Parameters
    ----------
    output_model: JWST data model
        flat-fielded input science data model, modified in-place

    flat_model: JWST data model
        data model containing flat-field

    Returns
    -------
    None
    """

    log.debug("Flat field correction for multi-slit mode.")

    any_updated = False # will set True if any flats applied

    # Check for a single image in the flat model
    if len(flat_model.slits) == 1 and flat_model.slits[0].name == None:
        # Populate slit attributes that will be needed later
        flat_model.slits[0].name = flat_model.meta.subarray.name
        flat_model.slits[0].xstart = flat_model.meta.subarray.xstart
        flat_model.slits[0].ystart = flat_model.meta.subarray.ystart
        flat_model.slits[0].xsize = flat_model.meta.subarray.xsize
        flat_model.slits[0].ysize = flat_model.meta.subarray.ysize

    # Apply flat to simple ImageModels
    if isinstance(output_model, datamodels.ImageModel):
        flat = get_flat(output_model, flat_model)
        if flat != None:
            apply_flat_field(output_model, flat)
            any_updated = True

    # Apply flat to MultiSlits
    elif isinstance(output_model, datamodels.MultiSlitModel):

        # Retrieve and apply flat to each slit contained in the input
        for slit in output_model.slits:
            log.info('Retrieving flat for slit %s' % (slit.name))
            flat = get_flat(slit, flat_model)
            if flat != None:
                apply_flat_field(slit, flat)
                any_updated = True

    # Apply flat to multiple-integration dataset
    elif isinstance(output_model, datamodels.CubeModel):
        sci_data = output_model.data
        nints = sci_data.shape[0]   # number of integrations

        #  Loop over integrations to flat field each
        for integ in range(nints):
            data_slice = sci_data[integ, :, :]

            # Create model for integration's data, w/needed subarray metadata
            integ_model = datamodels.ImageModel(data_slice)
            integ_model.update(output_model)

            log.info('Retrieving flat for integration %s' % (integ))
            flat = get_flat(integ_model, flat_model)
            if flat != None:
                apply_flat_field(integ_model, flat)
                any_updated = True

    if any_updated:
        output_model.meta.cal_step.flat_field = 'COMPLETE'
    else:
        output_model.meta.cal_step.flat_field = 'SKIPPED'


def get_flat(slit, flat_model):
    """
    Short Summary
    -------------
    For a simple ImageModel, get the single flat. Otherwise (for the multislit
    model), get the flat having the same name as the specified slit.

    Parameters
    ----------
    slit: JWST data model
        output data model

    flat_model: JWST data model
        data model containing flat-field

    Returns
    -------
    flat: 2D array
        flat field array for output model
    """

    if len(flat_model.slits) == 1:
        return flat_model.slits[0]
    else:
        found_one = False
        for flat in flat_model.slits:
            if flat.name == slit.name:
                log.info('Found matching flat %s' % (flat.name))
                found_one = True
                return flat

        if not found_one:
            log.error("Couldn't find matching flat-field for slit %s" % (slit.name))
            return None


def apply_flat_field(science, flat):
    """
    Short Summary
    -------------
    Flat fields the data and error arrays, and updates data quality array
    based on bad pixels in flat field arrays. Applies portion of flat field
    corresponding to science image subarray.

    Parameters
    ----------
    science: JWST data model
        input science data model

    flat: JWST data model
        flat field data model

    Returns
    -------
    None
    """

    # If the input science data model is a subarray, extract the same
    # subarray from the flatfield model
    if ref_matches_sci(flat, science):
        flat_data = flat.data
        flat_dq = flat.dq
    else:
        log.info("Extracting matching subarray from flat")
        flat_data = get_subarray(flat.data, science)
        flat_dq = get_subarray(flat.dq, science)

    # For pixels whose flat is either NaN or NO_FLAT_FIELD, update their DQ to
    # indicate that no flat is applied to those pixels
    flat_dq[np.isnan(flat_data)] = np.bitwise_or(flat_dq[np.isnan(flat_data)],
                                                 dqflags.pixel['NO_FLAT_FIELD'])

    # Replace NaN's in flat with 1's
    flat_data[np.isnan(flat_data)] = 1.0

    # Reset flat values of pixels having DQ values containing NO_FLAT_FIELD
    # to 1.0, so that no flat fielding correction is made
    wh_dq = np.bitwise_and(flat_dq, dqflags.pixel['NO_FLAT_FIELD'])
    flat_data[wh_dq == dqflags.pixel['NO_FLAT_FIELD']] = 1.0

    # Flatten data and error arrays
    science.data /= flat_data
    science.err /= flat_data

    # Combine the science and flat DQ arrays
    science.dq = np.bitwise_or(science.dq, flat_dq)


def ref_matches_sci(ref_model, sci_model):
    """
    Short Summary
    -------------
    Check if the science model has the same subarray parameters as the
    reference model.

    Parameters
    ----------
    ref_model: JWST data model
        data model containing flat-field

    sci_model: JWST data model
        input science data model to be flat-fielded

    Returns
    -------
    True if the science model has the same subarray parameters as the
    reference model, False otherwise.

    """
    # Get the science model subarray parameters
    try:
        xstart = sci_model.xstart
        xsize = sci_model.xsize
        ystart = sci_model.ystart
        ysize = sci_model.ysize
    except:
        xstart = sci_model.meta.subarray.xstart
        xsize = sci_model.meta.subarray.xsize
        ystart = sci_model.meta.subarray.ystart
        ysize = sci_model.meta.subarray.ysize

    log.debug(' sci xstart=%d, xsize=%d', xstart, xsize)
    log.debug(' sci ystart=%d, ysize=%d', ystart, ysize)
    log.debug(' ref xstart=%d, xsize=%d', ref_model.xstart, ref_model.xsize)
    log.debug(' ref ystart=%d, ysize=%d', ref_model.ystart, ref_model.ysize)

    # See if they match the reference model subarray parameters
    if (ref_model.xstart == xstart and ref_model.xsize == xsize and
        ref_model.ystart == ystart and ref_model.ysize == ysize):
        return True
    else:
        return False


def get_subarray(input_array, sci_model):
    """
    Short Summary
    -------------
    Return the slice from the input array using the subarray parameters of the
    science model.

    Parameters
    ----------
    input_array: 2D numpy array
        input array from which a subarray is extracted

    sci_model: JWST data model
        input science data model

    Returns
    -------
    A slice: 2D numpy array
        slice from the input array
    """

    # Get the science model subarray parameters
    try:
        xstart = sci_model.xstart
        xsize = sci_model.xsize
        ystart = sci_model.ystart
        ysize = sci_model.ysize
    except:
        xstart = sci_model.meta.subarray.xstart
        xsize = sci_model.meta.subarray.xsize
        ystart = sci_model.meta.subarray.ystart
        ysize = sci_model.meta.subarray.ysize

    # Compute the slicing indexes
    xstart = xstart - 1
    xstop = xstart + xsize
    ystart = ystart - 1
    ystop = ystart + ysize

    # Return the slice from the input array
    return input_array[ystart:ystop, xstart:xstop]

#
# The following functions are for multi-object spectroscopy flat fielding.
#
def do_MOS_flat_field(output_model, flat_model, flat_suffix=None):
    """
    Short Summary
    -------------
    Apply flat-fielding for MOS data, updating the output model.

    Parameters
    ----------
    output_model: JWST data model
        Science data model, modified (flat fielded) in-place.

    flat_model: JWST CubeFlatModel object
        Data model containing arrays for determining the flat field
        for each slit in `output_model`.

    flat_suffix: str or None
        Filename suffix for optional output file to save the interpolated
        flat field images.

    Returns
    -------
    MultiSlitModel or None
        If not None, the value will be the interpolated flat fields.

    """
    log.debug("Flat field correction for MOS data.")

    # All of these are 3-D arrays.
    flat_value = flat_model.data        # values of the flat field ...
    flat_wavelength = flat_model.wavelength     # ... at these wavelengths
    flat_data_qual = flat_model.dq

    # Are the wavelengths in the reference file increasing or decreasing?
    direction = find_dir(flat_wavelength)
    for_info = ["decreasing", "", "increasing"]
    log.info("Wavelengths in reference file are %s." % for_info[direction + 1])

    # Find pixels in the flat field that are NaN, and flag them as bad
    # in the data quality array.
    nan_bad = np.isnan(flat_value)
    flat_data_qual[nan_bad] |= dqflags.pixel['NO_FLAT_FIELD']
    del nan_bad

    # Reset flat values of pixels having DQ values containing NO_FLAT_FIELD
    # to 1.0, so that no flat field correction will be made.  This will
    # include all the pixels that were NaN in the flat field.
    bad = np.where(np.bitwise_and(flat_data_qual,
                                  dqflags.pixel['NO_FLAT_FIELD']))
    flat_value[bad] = 1.0
    del bad

    # Find pixels in the wavelength array that are zero, and flag them as
    # bad in the data quality array.  This is done after the check on the
    # NO_FLAT_FIELD flag in the section above, so as to avoid unnecessarily
    # modifying the flat field.
    zero_wl = (flat_wavelength == 0.)
    flat_data_qual[zero_wl] |= dqflags.pixel['NO_FLAT_FIELD']

    # Create an output model for the interpolated flat fields.
    if flat_suffix is not None:
        interpolated_flats = datamodels.MultiSlitModel()
        interpolated_flats.update(output_model, only="PRIMARY")
    else:
        interpolated_flats = None

    for (k, slit) in enumerate(output_model.slits):
        # Get the wavelength of each pixel in the extracted slit data.
        ysize, xsize = slit.data.shape
        grid = np.indices((ysize, xsize), dtype=np.float64)
        # The arguments are the X and Y pixel coordinates (in that order).
        (ra, dec, wl) = slit.meta.wcs(grid[1], grid[0])
        del ra, dec, grid

        xstart = slit.xstart - 1
        ystart = slit.ystart - 1
        xstop = xstart + xsize
        ystop = ystart + ysize

        # Extract the relevant slices of the reference file arrays.
        flat_val = flat_value[:, ystart:ystop, xstart:xstop]
        flat_wl = flat_wavelength[:, ystart:ystop, xstart:xstop]
        flat_dq = flat_data_qual[:, ystart:ystop, xstart:xstop]

        (flat_2d, flat_dq_2d) = interpolate_flat(flat_val, flat_dq,
                                                 flat_wl, wl, direction)
        if flat_suffix is not None:
            # Save flat_2d and flat_dq_2d for an output file.
            new_flat = datamodels.ImageModel(data=flat_2d, dq=flat_dq_2d)
            interpolated_flats.slits.append(new_flat.copy())
            interpolated_flats.slits[k].err[...] = 1.   # xxx not realistic
            nslit = len(interpolated_flats.slits) - 1
            interpolated_flats.slits[k].name = slit.name
            interpolated_flats.slits[k].xstart = slit.xstart
            interpolated_flats.slits[k].xsize = slit.xsize
            interpolated_flats.slits[k].ystart = slit.ystart
            interpolated_flats.slits[k].ysize = slit.ysize
            # Copy the WCS info from output (same as input).
            interpolated_flats.slits[k].meta.wcs = \
                  output_model.slits[k].meta.wcs

        slit.data /= flat_2d
        slit.err /= flat_2d
        slit.dq |= flat_dq_2d

        any_updated = True

    if any_updated:
        output_model.meta.cal_step.flat_field = 'COMPLETE'
    else:
        output_model.meta.cal_step.flat_field = 'SKIPPED'

    return interpolated_flats

def interpolate_flat(flat_val, flat_dq, flat_wl, wl, direction):
    """
    Short Summary
    -------------
    Interpolate within the 3-D flat field to get a 2-D flat.

    Parameters
    ----------
    flat_val: ndarray, 3-D
        This is slice [:, ystart:ystop, xstart:xstop] of the flat field
        reference image.  This slice covers the spatial extent of the
        extracted 2-D spectrum and includes all of the wavelength axis
        (the first axis) of the reference image.

    flat_dq: ndarray, 3-D
        This is slice [:, ystart:ystop, xstart:xstop] of the flat field
        data quality image.

    flat_wl: ndarray, 3-D
        The 3-D wavelength array for the flat field.  flat_wl[k, j, i] is
        the kth wavelength (at row j and column i) at which the flat field
        (in flat_val) was measured.  The wavelength wl[j, i] of the
        observed data at that pixel should lie within the range of values
        indexed by k.

    wl: ndarray, 2-D
        The wavelength at each pixel of the 2-D extracted spectrum.

    direction: int, +1 or -1
        +1 indicates that the wavelengths in flat_wl increase from one
        plane to the next, while -1 indicates that the wavelengths
        decrease.

    Returns
    -------
    tuple, two 2-D ndarrays
    flat_2d: ndarray, 2-D, float
        The flat field, interpolated over wavelength, same shape as `wl`.
        Divide the 2-D extracted spectrum by this array to correct for
        flat-field variations.
    flat_dq_2d: ndarray, 2-D, int
        The data quality array corresponding to flat_2d.
    """

    (nz, ysize, xsize) = flat_val.shape
    if nz == 1:
        return (flat_val[0], flat_dq[0])

    grid = np.indices((ysize, xsize), dtype=np.intp)
    ixpixel = grid[1]
    iypixel = grid[0]

    # The number of wavelengths at which the flat field was measured can
    # be as many as flat_wl.shape[0], but the actual number may be smaller;
    # call it num_k.  num_k can vary from pixel to pixel within flat_wl.
    temp = np.bitwise_and(flat_dq, dqflags.pixel['NO_FLAT_FIELD'])
    not_flagged = np.logical_not(temp.astype(np.bool))
    num_k = not_flagged.sum(axis=0, dtype=np.intp)
    max_num_k = num_k.max()
    if max_num_k <= 0:
        log.warning("Every pixel in the current slit is flagged with"
                    " NO_FLAT_FIELD.")
        return (flat_val[0], flat_dq[0])

    # Maximum index for the first axis; this can vary from pixel to pixel.
    k_max = num_k - 1
    k_max = np.where(k_max < 0, 0, k_max)       # k_max must not be negative
    k_max_max = k_max.max()                     # for the range of a loop

    # Arrays for the flat field, wavelengths, and data quality array,
    # but compressed (see below) along the first axis to remove all values
    # that are flagged with NO_FLAT_FIELD in the data quality array.
    compr_val = np.ones((max_num_k, ysize, xsize), dtype=flat_val.dtype)
    compr_wl = np.zeros((max_num_k, ysize, xsize), dtype=flat_wl.dtype)
    compr_dq = np.zeros((max_num_k, ysize, xsize), dtype=flat_dq.dtype)

    for j in range(ysize):
        for i in range(xsize):
            compr_val[0:num_k[j, i], j, i] = np.compress(not_flagged[:, j, i],
                                                        flat_val[:, j, i],
                                                        axis=0)
            compr_wl[0:num_k[j, i], j, i] = np.compress(not_flagged[:, j, i],
                                                       flat_wl[:, j, i],
                                                       axis=0)
            compr_dq[0:num_k[j, i], j, i] = np.compress(not_flagged[:, j, i],
                                                       flat_dq[:, j, i],
                                                       axis=0)

    # If there's no good data at a pixel, flag that pixel with NO_FLAT_FIELD.
    # The input probably was already flagged, but if num_k is zero, that
    # flag would have been lost by calling np.compress() just above.
    compr_dq[0, :, :] = np.where(num_k == 0,
                                 dqflags.pixel['NO_FLAT_FIELD'],
                                 compr_dq[0, :, :])

    # The initial value of -1 is a flag to indicate that elements have not
    # been assigned valid values yet.
    k = np.zeros(wl.shape, dtype=np.intp) - 1

    # Truncate the index for wavelengths that are outside the range of
    # flat_wl.  Near the end of this function we'll flag these pixels
    # with NO_FLAT_FIELD and set flat_2d to 1, but for now the indices
    # need to be assigned harmless values to avoid indexing out of bounds.
    # Note that k_max is a 2-D array, one element for each image pixel.
    #   Why set the upper limit of k to k_max - 1 instead of just k_max?
    #   Because we interpolate using elements k and k + 1.
    if direction > 0:
        k[:, :] = np.where(wl <= compr_wl[0, :, :], 0, k)
        k[:, :] = np.where(wl >= compr_wl[k_max, iypixel, ixpixel],
                           k_max - 1, k)
    else:
        k[:, :] = np.where(wl >= compr_wl[0, :, :], 0, k)
        k[:, :] = np.where(wl <= compr_wl[k_max, iypixel, ixpixel],
                           k_max - 1, k)

    # Look for the correct interval for linear interpolation.
    for k_test in range(k_max_max):
        if direction > 0:
            test1 = np.logical_and(wl >= compr_wl[k_test, :, :],
                                   wl < compr_wl[k_test + 1, :, :])
        else:
            test1 = np.logical_and(wl <= compr_wl[k_test, :, :],
                                   wl > compr_wl[k_test + 1, :, :])
        # If an element of k is not -1, it has already been assigned, and
        # I don't want to clobber it.
        test2 = np.logical_and(k == -1, test1)
        k[:, :] = np.where(test2, k_test, k)
        if np.all(k >= 0):
            break
    del test1, test2

    # At this point, all elements of k should have been assigned a value,
    # but check to be sure.
    if np.any(k == -1):
        n_neg = (k == -1).sum(dtype=np.int64)
        log.error("Internal error:  %d elements of k are -1" % n_neg)
        k[:, :] = np.where(k == -1, 0, k)
        del n_neg

    # k[j,i] is the element in compr_wl such that
    # compr_wl[k[j,i], j, i] <= wl[j,i] <= compr_wl[k[j,i]+1, j, i]
    # or >=, depending on `direction`.

    # Use linear interpolation within the 3-D flat field to get a 2-D
    # flat field.
    denom = compr_wl[k + 1, iypixel, ixpixel] - compr_wl[k, iypixel, ixpixel]
    zero_denom = (denom == 0.)
    denom = np.where(zero_denom, 1., denom)
    q = np.where(zero_denom, 0.,
                 (wl - compr_wl[k, iypixel, ixpixel]) / denom)
    p = 1. - q
    flat_2d = p * compr_val[k, iypixel, ixpixel] + \
              q * compr_val[k + 1, iypixel, ixpixel]
    flat_dq_2d = np.where(q == 0.,
                          compr_dq[k, iypixel, ixpixel],
                          np.bitwise_or(compr_dq[k, iypixel, ixpixel],
                                        compr_dq[k + 1, iypixel, ixpixel]))
    del p, q, denom

    # If the wavelength is out of range, set the 2-D flat field to 1,
    # and set a flag in the output 2-D DQ array.
    if direction > 0:
        indx = np.where(wl < compr_wl[0, :, :])
        flat_2d[indx] = 1.
        flat_dq_2d[indx] |= dqflags.pixel['NO_FLAT_FIELD']
        indx = np.where(wl > compr_wl[k_max, iypixel, ixpixel])
        flat_2d[indx] = 1.
        flat_dq_2d[indx] |= dqflags.pixel['NO_FLAT_FIELD']
    else:
        indx = np.where(wl > compr_wl[0, :, :])
        flat_2d[indx] = 1.
        flat_dq_2d[indx] |= dqflags.pixel['NO_FLAT_FIELD']
        indx = np.where(wl < compr_wl[k_max, iypixel, ixpixel])
        flat_2d[indx] = 1.
        flat_dq_2d[indx] |= dqflags.pixel['NO_FLAT_FIELD']

    return (flat_2d.astype(flat_val.dtype), flat_dq_2d)

def find_dir(flat_wavelength):
    """Check whether the values are increasing or decreasing

    Parameters
    ----------
    flat_wavelength, ndarray
        The 3-D array of wavelengths as read from the flat field
        reference file.  Whether these are increasing or decreasing will
        be determined by comparing the first and second planes.

    Returns
    -------
    int, -1 or +1
        The returned value will be +1 if the values in flat_wavelength
        increase as the first index increases (or if there is only one
        plane, or no non-zero values in common to both the first and
        second planes); otherwise, the returned value will be -1.
    """

    if len(flat_wavelength) < 2:
        return 1

    # The comparison would not be meaningful if some elements were zero
    # in one plane but non-zero in the other.
    nonzero0 = (flat_wavelength[0, ...] > 0.)
    nonzero1 = (flat_wavelength[1, ...] > 0.)
    nonzero = np.logical_and(nonzero0, nonzero1)

    # diff is 1-D
    diff = flat_wavelength[1, nonzero] - flat_wavelength[0, nonzero]
    nelem = diff.shape[0]
    if nelem < 1:
        return 1

    diff.sort()
    if diff[nelem // 2] >= 0:
        dir = 1
    else:
        dir = -1

    return dir
