#
#  Module for IPC correction.
#

from collections import namedtuple
import logging
import numpy as np

from ..lib import pipe_utils
from . import x_irs2

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

NumRefPixels = namedtuple("NumRefPixels",
                          ["bottom_rows", "top_rows",
                           "left_columns", "right_columns"])


def do_correction(input_model, ipc_model):
    """Execute all tasks for IPC correction

    Parameters
    ----------
    input_model : data model object
        Science data to be corrected.

    ipc_model : ipc model object
        Deconvolution kernel, either a 2-D or 4-D image in the first
        extension.

    Returns
    -------
    output_model : data model object
        IPC-corrected science data.
    """

    # Save some data params for easy use later
    sci_nints = input_model.data.shape[0]
    sci_ngroups = input_model.data.shape[1]
    sci_nframes = input_model.meta.exposure.nframes
    sci_groupgap = input_model.meta.exposure.groupgap

    log.debug('IPC corr using nints=%d, ngroups=%d, nframes=%d, groupgap=%d' %
              (sci_nints, sci_ngroups, sci_nframes, sci_groupgap))

    # Apply the correction.
    output_model = ipc_correction(input_model, ipc_model)

    return output_model


def ipc_correction(input_model, ipc_model):
    """Apply the IPC correction to the science arrays.

    Parameters
    ----------
    input_model : data model object
        The input science data.

    ipc_model : IPCModel object
        The IPC kernel.  The input is corrected for IPC by convolving
        with this 2-D or 4-D array.

    Returns
    -------
    output : data model object
        IPC-corrected science data.
    """

    log.debug("ipc_correction: nints=%d, ngroups=%d, size=%d,%d",
              input_model.meta.exposure.nints,
              input_model.meta.exposure.ngroups,
              input_model.data.shape[-1],
              input_model.data.shape[-2])

    # Create output as a copy of the input science data model.
    output = input_model.copy()

    # Was IRS2 readout used?
    is_irs2_format = pipe_utils.is_irs2(input_model)
    if is_irs2_format:
        irs2_mask = x_irs2.make_mask(input_model)

    detector = input_model.meta.instrument.detector

    # The number of reference pixels along the bottom edge, top edge,
    # left edge, and right edge.
    nref = get_num_ref_pixels(input_model)

    # Get the data for the IPC kernel.  This can be a slice, if input_model
    # is a subarray.
    kernel = get_ipc_slice(input_model, ipc_model)

    log.debug("substrt1 = %d, subsize1 = %d, substrt2 = %d, subsize2 = %d" %
              (input_model.meta.subarray.xstart, input_model.meta.subarray.xsize,
               input_model.meta.subarray.ystart, input_model.meta.subarray.ysize))
    log.debug('Number of reference pixels: bottom, top, left, right ='
              ' %d, %d, %d, %d' %
              (nref.bottom_rows, nref.top_rows,
               nref.left_columns, nref.right_columns))
    log.debug("Shape of ipc image = %s" % repr(ipc_model.data.shape))

    # Loop over all integrations and groups in input science data.
    for i in range(input_model.data.shape[0]):                  # integrations
        for j in range(input_model.data.shape[1]):              # groups
            # Convolve the current group in-place with the IPC kernel.
            if is_irs2_format:
                # Extract normal data from input IRS2-format data.
                temp = x_irs2.from_irs2(output.data[i, j, :, :], irs2_mask,
                                        detector)
                ipc_convolve(temp, kernel, nref)
                # Insert normal data back into original, IRS2-format data.
                x_irs2.to_irs2(output.data[i, j, :, :], temp, irs2_mask,
                               detector)
            else:
                ipc_convolve(output.data[i, j], kernel, nref)

    return output


def get_num_ref_pixels(input_model):
    """Get the number of reference pixel rows and columns.

    Parameters
    ----------
    input_model : data model object
        The input science data.

    Returns
    -------
    nref : namedtuple
        bottom_rows : int
            The number of reference rows at the bottom of the image.
        top_rows : int
            The number of reference rows at the top of the image.
        left_columns : int
            The number of reference columns at the left edge.
        right_columns : int
            The number of reference columns at the right edge.
    """

    xstart = input_model.meta.subarray.xstart - 1       # zero indexed
    xsize = input_model.meta.subarray.xsize
    if input_model.meta.instrument.name == 'MIRI':
        nref = NumRefPixels(bottom_rows=0,
                            top_rows=0,
                            left_columns=max(0, 4 - xstart),
                            right_columns=max(0, xstart + xsize - 1028))
    else:
        ystart = input_model.meta.subarray.ystart - 1   # zero indexed
        ysize = input_model.meta.subarray.ysize
        nref = NumRefPixels(bottom_rows=max(0, 4 - ystart),
                            top_rows=max(0, ystart + ysize - 2044),
                            left_columns=max(0, 4 - xstart),
                            right_columns=max(0, xstart + xsize - 2044))

    return nref


def get_ipc_slice(input_model, ipc_model):
    """Extract a slice from IPC kernel corresponding to science data.

    Parameters
    ----------
    input_model : data model object
        The input science data.

    ipc_model : data model object
        The IPC kernel model.

    Returns
    -------
    kernel : ndarray, either 2-D or 4-D
        The data array for the IPC kernel.  If the IPC kernel is 4-D and
        the science data array is a subarray, `kernel` will be a slice of
        the reference image; otherwise, this will be the full image.
    """

    if len(ipc_model.data.shape) == 2:
        return ipc_model.data

    # Convert xstart and ystart from one indexing to zero indexing.
    xstart = input_model.meta.subarray.xstart - 1
    ystart = input_model.meta.subarray.ystart - 1
    xsize = input_model.meta.subarray.xsize
    ysize = input_model.meta.subarray.ysize
    if input_model.meta.instrument.name == 'MIRI':
        is_subarray = (xsize < 1032 or ysize < 1024)
    else:
        is_subarray = (xsize < 2048 or ysize < 2048)

    if is_subarray:
        return ipc_model.data[:, :, ystart:ystart + ysize, xstart:xstart + xsize]
    else:
        return ipc_model.data


def ipc_convolve(output_data, kernel, nref):
    """Convolve the science data with the IPC kernel.

    Parameters
    ----------
    output_data : ndarray, 2-D
        A copy of the input science data for one group; this will be
        modified in-place.

    kernel : ndarray, 2-D or 4-D
        The IPC kernel; the input is corrected for IPC by convolving with
        this array.  If it is 4-D the last two dimensions will be a slice
        that matches the last two dimensions of `output_data`.

    nref : namedtuple
        bottom_rows : int
            The number of reference rows at the bottom of the image.
        top_rows : int
            The number of reference rows at the top of the image.
        left_columns : int
            The number of reference columns at the left edge.
        right_columns : int
            The number of reference columns at the right edge.
    """

    bottom_rows = nref.bottom_rows
    top_rows = nref.top_rows
    left_columns = nref.left_columns
    right_columns = nref.right_columns

    kshape = kernel.shape

    # This is the shape of the entire image, which may include reference
    # pixels.
    shape = output_data.shape

    # These axis lengths exclude reference pixels, if there are any.
    ny = shape[0] - (bottom_rows + top_rows)
    nx = shape[1] - (left_columns + right_columns)

    # The temporary array temp is larger than the science part of
    # output_data by a border (set to zero) that's about half of the
    # kernel size, so the convolution can be done without checking for
    # out of bounds.
    # b_b, t_b, l_b, and r_b are the widths of the borders on the
    # bottom, top, left, and right, respectively.
    b_b = kshape[0] // 2
    t_b = kshape[0] - b_b - 1
    l_b = kshape[1] // 2
    r_b = kshape[1] - l_b - 1
    tny = ny + b_b + t_b
    yoff = bottom_rows                      # offset in output_data
    tnx = nx + l_b + r_b
    xoff = left_columns                     # offset in output_data

    # Note that when we accumulate sums to output_data below, we will
    # always use the same slice:  output_data[yoff:yoff+ny, xoff:xoff+nx].

    # Copy the science portion (not the reference pixels) of output_data
    # to this temporary array, then make subsequent changes in-place to
    # output_data.
    temp = np.zeros((tny, tnx), dtype=output_data.dtype)
    temp[b_b:b_b + ny, l_b:l_b + nx] = \
        output_data[yoff:yoff + ny, xoff:xoff + nx].copy()

    # After setting this slice to zero, we'll incrementally add to it.
    output_data[yoff:yoff + ny, xoff:xoff + nx] = 0.

    if len(kshape) == 2:
        # 2-D IPC kernel.  Loop over pixels of the deconvolution kernel.
        # In this section, `part` has the same shape as `temp`.
        middle_j = kshape[0] // 2
        middle_i = kshape[1] // 2
        for j in range(kshape[0]):
            jstart = kshape[0] - j - 1
            for i in range(kshape[1]):
                if i == middle_i and j == middle_j:
                    continue                # the middle pixel is done last
                part = kernel[j, i] * temp
                istart = kshape[1] - i - 1
                output_data[yoff:yoff + ny, xoff:xoff + nx] += \
                    part[jstart:jstart + ny, istart:istart + nx]
        # The middle pixel of the IPC kernel is expected to be the largest,
        # so add that last.
        part = kernel[middle_j, middle_i] * temp
        output_data[yoff:yoff + ny, xoff:xoff + nx] += \
            part[middle_j:middle_j + ny, middle_i:middle_i + nx]

    else:
        # 4-D IPC kernel.  Extract a subset of the kernel:  all of the
        # first two axes, but only the portion of the last two axes
        # corresponding to the science data (i.e. possibly a subarray,
        # and certainly excluding reference pixels).
        k_temp = np.zeros((kshape[0], kshape[1], tny, tnx),
                          dtype=kernel.dtype)
        k_temp[:, :, b_b:b_b + ny, l_b:l_b + nx] = \
            kernel[:, :, yoff:yoff + ny, xoff:xoff + nx]

        # In this section, `part` has shape (ny, nx), which is smaller
        # than `temp`.
        middle_j = kshape[0] // 2
        middle_i = kshape[1] // 2
        for j in range(kshape[0]):
            jstart = kshape[0] - j - 1
            for i in range(kshape[1]):
                if i == middle_i and j == middle_j:
                    continue                # the middle pixel is done last
                istart = kshape[1] - i - 1
                # The slice of k_temp includes different pixels for the
                # first or second axes within each loop, but the same slice
                # for the last two axes.
                # The slice of temp (a copy of the science data) includes
                # a different offset for each loop.
                part = k_temp[j, i, b_b:b_b + ny, l_b:l_b + nx] * \
                    temp[jstart:jstart + ny, istart:istart + nx]
                output_data[yoff:yoff + ny, xoff:xoff + nx] += part
        # Add the product for the middle pixel last.
        part = k_temp[middle_j, middle_i, b_b:b_b + ny, l_b:l_b + nx] * \
            temp[middle_j:middle_j + ny, middle_i:middle_i + nx]
        output_data[yoff:yoff + ny, xoff:xoff + nx] += part
