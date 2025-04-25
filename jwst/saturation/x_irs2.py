from collections import namedtuple

import numpy as np
from ..lib import pipe_utils

""" This is the interface:
    mask = make_mask(input_model)
        Create a mask for extracting normal pixels; used by from_irs2 and
        to_irs2.
        n and r can be specified as keyword arguments to override the
        default values.  This option is primarily useful if `input_model`
        is a numpy.ndarray rather than a jwst.datamodels object.
    shape = normal_shape(input_model, n=n, r=r)
        The shape of the data array when excluding interleaved reference
        pixels.
        n and r can be specified as keyword arguments.
    normal_data = from_irs2(irs2_data, mask, detector)
        Extract the normal pixels from data in IRS2 format.
    to_irs2(irs2_data, normal_data, mask, detector)
        Insert an array of normal pixels back into data in IRS2 format.

    Note that `input_model` may be either a jwst.datamodels object or a
    numpy.ndarray (though in the latter case the parameters will be
    assigned default values, unless specified explicitly).
"""


ReadoutParam = namedtuple("ReadoutParam", ["refout", "n", "r"])


def _get_irs2_parameters(input_model, n=None, r=None):
    """
    Get the parameters describing IRS2 readout format.

    Parameters
    ----------
    input_model : DataModel or numpy.ndarray
        Model from which we retrieve the width of the reference
        output and the values of NRS_NORM and NRS_REF. If the input_model
        is a ndarray the parameters will be assigned default values.

    n : int or None
        The number of "normal" (as opposed to reference) pixels read
        out before jumping to the reference pixel region.

    r : int or None
        The number of reference pixels read out before jumping back to
        the normal pixel region.

    Returns
    -------
    param : ReadoutParam object
        ReadoutParam objects contains the following information:
        param.refout: int
            The length (in the last image axis) of the reference output
            section.  The reference output is assumed to be on the left
            side of the IRS2-format image.

        param.n : int
            The number of "normal" (as opposed to reference) pixels read
            out before jumping to the reference pixel region.

        param.r : int
            The number of reference pixels read out before jumping back to
            the normal pixel region.
    """
    try:
        # Try to get keyword values
        n_norm = input_model.meta.exposure.nrs_normal
        n_ref = input_model.meta.exposure.nrs_reference
    except AttributeError:
        # If keywords are missing, use default values
        n_norm = 16
        n_ref = 4

    # Check for user-supplied values
    if n is not None:
        n_norm = n
    if r is not None:
        n_ref = r

    param = ReadoutParam(refout=(512 + 512 // n_norm * n_ref), n=n_norm, r=n_ref)

    return param


def normal_shape(input_model, n=None, r=None, detector=None):
    """
    Determine the shape of the 'normal' pixel data.

    Parameters
    ----------
    input_model : DataModel
        This is used to get the shape of the input data.

    n : int or None
        The number of "normal" (as opposed to reference) pixels read
        out before jumping to the reference pixel region.

    r : int or None
        The number of reference pixels read out before jumping back to
        the normal pixel region.

    detector : str
       Detector of data. Valid values are None, NRS1, or NRS2. Other detector
       values will result in a run time error.

    Returns
    -------
    data_shape : 2 D array
        The shape of the data array when excluding interleaved reference
        pixels.
    """
    if isinstance(input_model, np.ndarray):
        shape = input_model.shape
    else:
        shape = input_model.data.shape
        if detector is None:
            detector = input_model.meta.instrument.detector

    if not pipe_utils.is_irs2(input_model):  # not IRS2 format
        return shape

    param = _get_irs2_parameters(input_model, n=n, r=r)

    if detector is None:
        irs2_nx = shape[-1]
    elif detector == "NRS1" or detector == "NRS2":
        irs2_nx = shape[-2]
    else:
        raise RuntimeError(f"Detector {detector} is not supported for IRS2 data.")

    k = (irs2_nx - param.refout) // (param.n + param.r)
    n_output = (irs2_nx - param.refout) - k * param.r

    if detector is None:
        data_shape = shape[0:-1] + (n_output,)
    elif detector == "NRS1" or detector == "NRS2":
        data_shape = shape[0:-2] + (n_output, shape[-1])

    return data_shape


def make_mask(input_model, n=None, r=None):
    """
    Create a mask to extract "normal" pixels.

    Parameters
    ----------
    input_model : the model for the input data, or a numpy.ndarray
        This is used for getting the IRS2 parameters and the length of
        the X image axis.

    n : int or None
        The number of "normal" (as opposed to reference) pixels read
        out before jumping to the reference pixel region.

    r : int or None
        The number of reference pixels read out before jumping back to
        the normal pixel region.

    Returns
    -------
    irs2_mask : 1-D bool array
        Boolean index mask with length equal to the last axis of
        the data shape.
    """
    param = _get_irs2_parameters(input_model, n=n, r=r)
    refout = param.refout
    n_norm = param.n
    n_ref = param.r

    if isinstance(input_model, np.ndarray):
        shape = input_model.shape
    else:
        shape = input_model.data.shape
    # The input may be flipped and/or rotated from detector orientation.
    irs2_nx = max(shape[-1], shape[-2])

    # Number of (n + r) per output, assuming 4 amplifier outputs.
    k = (irs2_nx - refout) // 4 // (n_norm + n_ref)
    # Number of normal pixels per amplifier output.
    n_output = (irs2_nx - refout) // 4 - k * n_ref

    irs2_mask = np.ones(irs2_nx, dtype=bool)
    irs2_mask[0:refout] = False

    # Check that the locations of interspersed reference pixels is
    # the same regardless of readout direction.
    if n_output // n_norm * n_norm == n_output:
        # The interspersed reference pixels are in the same locations
        # regardless of readout direction.
        for i in range(refout + n_norm // 2, irs2_nx + 1, n_norm + n_ref):
            irs2_mask[i : i + n_ref] = False
    else:
        # Set the flags for each readout direction separately.
        nelem = (irs2_nx - refout) // 4  # number of elements per output
        temp = np.ones(nelem, dtype=bool)
        for i in range(n_norm // 2, nelem + 1, n_norm + n_ref):
            temp[i : i + n_ref] = False
        j = refout
        irs2_mask[j : j + nelem] = temp.copy()
        j = refout + nelem
        irs2_mask[j + nelem - 1 : j - 1 : -1] = temp.copy()
        j = refout + 2 * nelem
        irs2_mask[j : j + nelem] = temp.copy()
        j = refout + 3 * nelem
        irs2_mask[j + nelem - 1 : j - 1 : -1] = temp.copy()

    return irs2_mask


def from_irs2(irs2_data, irs2_mask, detector=None):
    """
    Extract 'normal' pixel data from an IRS2 image.

    Parameters
    ----------
    irs2_data : ndarray
        Data in IRS2 format.  This can be a slice in the Y direction, but
        it should include the entire X (last) axis.

    irs2_mask : 1-D array, bool
        Boolean mask to extract the "normal" pixels.  This is a 1-D array
        with length equal to the size of the next-to-last axis (for data
        in DMS orientation) of `irs2_data`.

    detector : str or None
        For IRS2 data in DMS orientation, `detector` should be either
        "NRS1" or "NRS2"; NIRSpec is currently the only instrument
        supported in this module.  The mask will be applied to the rows,
        and for "NRS2" the mask will first be reversed.
        For IRS2 data in detector orientation, `detector` should be None
        (the default), and the mask will be applied to the columns.

    Returns
    -------
    norm_data : ndarray
       The normal pixel data (i.e. without embedded reference pixels).
    """
    if detector is None:
        # Select columns.
        norm_data = irs2_data[..., irs2_mask]
    elif detector == "NRS1":
        # Select rows.
        norm_data = irs2_data[..., irs2_mask, :]
    elif detector == "NRS2":
        # Reverse the direction of the mask, and select rows.
        temp_mask = irs2_mask[::-1]
        norm_data = irs2_data[..., temp_mask, :]
    else:
        raise RuntimeError(f"Detector {detector} is not supported for IRS2 data.")

    return norm_data


def to_irs2(irs2_data, norm_data, irs2_mask, detector=None):
    """
    Copy 'normal' pixel data into an IRS2 image.

    Parameters
    ----------
    irs2_data : ndarray
        Data in IRS2 format.  This will be modified in-place.

    norm_data : ndarray
        The normal data, for example previously extracted from irs2_data
        but then modified in some way.  This will be copied back into
        irs2_data in the correct locations, as specified by `irs2_mask`.

    irs2_mask : 1-D array, bool
        Boolean mask identifying the locations of the "normal" pixels
        within irs2_data.  The length is equal to the size of the
        next-to-last axis (for data in DMS orientation) of `irs2_data`.

    detector : str or None
        For IRS2 data in DMS orientation, `detector` should be either
        "NRS1" or "NRS2"; NIRSpec is currently the only instrument
        supported in this module.  The mask will be applied to the rows,
        and for "NRS2" the mask will first be reversed.
        For IRS2 data in detector orientation, `detector` should be None
        (the default), and the mask will be applied to the columns.
    """
    if detector is None:
        # Mask specifies columns.
        irs2_data[..., irs2_mask] = norm_data.copy()
    elif detector == "NRS1":
        # Mask specifies rows.
        irs2_data[..., irs2_mask, :] = norm_data.copy()
    elif detector == "NRS2":
        # Reverse the direction of the mask, and apply the reversed mask
        # to the rows.
        temp_mask = irs2_mask[::-1]
        irs2_data[..., temp_mask, :] = norm_data.copy()
    else:
        raise RuntimeError(f"Detector {detector} is not supported for IRS2 data.")
