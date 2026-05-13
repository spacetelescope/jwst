import logging

import numpy as np
from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)

__all__ = [
    "is_subarray",
    "ref_matches_sci",
    "get_subarray_model",
    "get_multistripe_subarray_model",
    "science_detector_frame_transform",
    "detector_science_frame_transform",
    "MatchRowError",
    "find_row",
]


def is_subarray(input_model):
    """
    Check if a data model comes from a subarray readout.

    If the data dimensions are less than 2048x2048 (or 1032x1024
    for MIRI), it is assumed to be a subarray.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Input data model to be checked.

    Returns
    -------
    status : bool
        `True` if the primary data array of the model is smaller than
        full-frame, `False` otherwise.
    """
    # Get the dimensions of the model's primary data array
    nrows = input_model.data.shape[-2]
    ncols = input_model.data.shape[-1]

    # Compare the dimensions to a full-frame
    if input_model.meta.instrument.name.upper() == "MIRI":
        if ncols < 1032 or nrows < 1024:
            log.debug("Input exposure is a subarray readout")
            return True
        else:
            return False
    else:
        if ncols < 2048 or nrows < 2048:
            log.debug("Input exposure is a subarray readout")
            return True
        else:
            return False


def ref_matches_sci(sci_model, ref_model):
    """
    Match science model to reference model.

    Check to see if the science model has the same subarray
    characteristics as the reference model. Also performs error
    checking on the subarray keywords in both the reference file
    and science data.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Science data model.

    ref_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Reference file data model.

    Returns
    -------
    status : bool
        `True` if the science model has the same subarray parameters
        as the reference model, `False` otherwise.
    """
    # Check first for superstripe data: it should never match the ref model.
    # Reference files for superstripe ramps are handled directly by
    # get_multistripe_subarray_model.
    nstripe = getattr(sci_model.meta.subarray, "num_superstripe", None)
    if nstripe is not None and nstripe > 0:
        return False

    # Get the reference file subarray parameters
    xstart_ref = ref_model.meta.subarray.xstart
    xsize_ref = ref_model.meta.subarray.xsize
    ystart_ref = ref_model.meta.subarray.ystart
    ysize_ref = ref_model.meta.subarray.ysize

    # Make sure the attributes were populated
    if (xstart_ref is None) or (xsize_ref is None) or (ystart_ref is None) or (ysize_ref is None):
        # If the ref file is full-frame, set the missing params to
        # default values
        if ref_model.meta.instrument.name.upper() == "MIRI":
            if ref_model.data.shape[-1] == 1032 and ref_model.data.shape[-2] == 1024:
                log.warning("Missing subarray corner/size keywords in reference file")
                log.warning("Setting them to full-frame default values")
                xstart_ref = 1
                ystart_ref = 1
                xsize_ref = 1032
                ysize_ref = 1024
                ref_model.meta.subarray.xstart = xstart_ref
                ref_model.meta.subarray.ystart = ystart_ref
                ref_model.meta.subarray.xsize = xsize_ref
                ref_model.meta.subarray.ysize = ysize_ref
            else:
                log.error("Missing subarray corner/size keywords in reference file")
                raise ValueError("Can't determine ref file subarray properties")
        else:
            if ref_model.data.shape[-1] == 2048 and ref_model.data.shape[-2] == 2048:
                log.warning("Missing subarray corner/size keywords in reference file")
                log.warning("Setting them to full-frame default values")
                xstart_ref = 1
                ystart_ref = 1
                xsize_ref = 2048
                ysize_ref = 2048
                ref_model.meta.subarray.xstart = xstart_ref
                ref_model.meta.subarray.ystart = ystart_ref
                ref_model.meta.subarray.xsize = xsize_ref
                ref_model.meta.subarray.ysize = ysize_ref
            else:
                log.error("Missing subarray corner/size keywords in reference file")
                raise ValueError("Can't determine ref file subarray properties")

    log.debug(
        "ref substrt1=%d, subsize1=%d, substrt2=%d, subsize2=%d",
        xstart_ref,
        xsize_ref,
        ystart_ref,
        ysize_ref,
    )

    # Make sure the starting corners are valid
    if xstart_ref < 1 or ystart_ref < 1:
        log.error("Reference file subarray corners are invalid")
        raise ValueError("Bad subarray corners")

    # Make sure the size attributes are consistent with data array
    if hasattr(ref_model, "coeffs"):
        xsize_data = ref_model.coeffs.shape[-1]
        ysize_data = ref_model.coeffs.shape[-2]
    else:
        xsize_data = ref_model.shape[-1]
        ysize_data = ref_model.shape[-2]

    # Check x dimensions
    if xsize_ref != xsize_data:
        log.warning("Reference file data array size doesn't match SUBSIZE1")
        log.warning("Using actual array size")
        xsize_ref = xsize_data
    # Check y dimensions
    if ysize_ref != ysize_data:
        # NIRSpec IRS2 is a special mode, where it's allowed to have
        # a mismatch in the y-size.
        if ysize_ref == 2048 and ysize_data == 3200:
            pass
        else:
            log.warning("Reference file data array size doesn't match SUBSIZE2")
            log.warning("Using actual array size")
            ysize_ref = ysize_data

    # Get the science model subarray parameters
    try:
        # If the science data model has already gone through
        # extract_2d, the metadata is in the top level of the model
        xstart_sci = sci_model.xstart
        xsize_sci = sci_model.xsize
        ystart_sci = sci_model.ystart
        ysize_sci = sci_model.ysize
    except AttributeError:
        # Otherwise the metadata is in the meta tree
        xstart_sci = sci_model.meta.subarray.xstart
        xsize_sci = sci_model.meta.subarray.xsize
        ystart_sci = sci_model.meta.subarray.ystart
        ysize_sci = sci_model.meta.subarray.ysize

    # Make sure the attributes were populated
    if (xstart_sci is None) or (xsize_sci is None) or (ystart_sci is None) or (ysize_sci is None):
        log.error("Missing subarray corner/size keywords in science file")
        raise ValueError("Can't determine science file subarray properties")
    else:
        log.debug(
            "sci substrt1=%d, subsize1=%d, substrt2=%d, subsize2=%d",
            xstart_sci,
            xsize_sci,
            ystart_sci,
            ysize_sci,
        )

    # Make sure the starting corners are valid
    if xstart_sci < 1 or ystart_sci < 1:
        log.error("Science file subarray corners are invalid")
        raise ValueError("Bad subarray corners")

    # Make sure the size attributes are consistent with data array.
    # Check x dimensions
    if xsize_sci != sci_model.shape[-1]:
        log.warning("Science file data array size doesn't match SUBSIZE1")
        log.warning("Using actual array size")
        xsize_sci = sci_model.shape[-1]
    # Check y dimensions
    if ysize_sci != sci_model.shape[-2]:
        # NIRSpec IRS2 is a special mode, where it's allowed to have
        # a mismatch in the y-size.
        if ysize_sci == 2048 and sci_model.shape[-2] == 3200:
            pass
        else:
            log.warning("Science file data array size doesn't match SUBSIZE2")
            log.warning("Using actual array size")
            ysize_sci = sci_model.shape[-2]

    # Finally, see if all of the science file subarray params match those
    # of the reference file
    if (
        xstart_ref == xstart_sci
        and xsize_ref == xsize_sci
        and ystart_ref == ystart_sci
        and ysize_ref == ysize_sci
    ):
        return True
    else:
        return False


def get_subarray_model(sci_model, ref_model):
    """
    Get subarray model.

    Create a subarray version of a reference file model that matches
    the subarray characteristics of a science data model. A new
    model is created that contains subarrays of all data arrays
    contained in the reference file model.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Science data model.

    ref_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Reference file data model.

    Returns
    -------
    sub_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Subarray version of the reference file model.
    """
    # If science data is in multistripe readout, use
    # multistripe-specific subarray reconstruction.
    if not isinstance(sci_model, datamodels.JwstDataModel):
        raise TypeError("Science model must be a JWST data model.")
    if not isinstance(ref_model, datamodels.JwstDataModel):
        raise TypeError("Reference model must be a JWST data model.")
    if getattr(sci_model.meta.subarray, "multistripe_reads1", None) is not None:
        return get_multistripe_subarray_model(sci_model, ref_model)

    # Get the science model subarray params
    xstart_sci = sci_model.meta.subarray.xstart
    xsize_sci = sci_model.meta.subarray.xsize
    ystart_sci = sci_model.meta.subarray.ystart
    ysize_sci = sci_model.meta.subarray.ysize

    # Get the reference model subarray params
    xstart_ref = ref_model.meta.subarray.xstart
    ystart_ref = ref_model.meta.subarray.ystart
    xsize_ref = ref_model.meta.subarray.xsize
    ysize_ref = ref_model.meta.subarray.ysize

    # Compute the slice indexes, in 0-indexed python frame
    xstart = xstart_sci - xstart_ref
    ystart = ystart_sci - ystart_ref
    xstop = xstart + xsize_sci
    ystop = ystart + ysize_sci
    log.debug("slice xstart=%d, xstop=%d, ystart=%d, ystop=%d", xstart, xstop, ystart, ystop)

    # Make sure that the slice limits are within the bounds of
    # the reference file data array
    if (
        xstart < 0
        or ystart < 0
        or xstop > ref_model.meta.subarray.xsize
        or ystop > ref_model.meta.subarray.ysize
    ):
        log.error(
            "Computed reference file slice indexes are incompatible with "
            "size of reference data array"
        )
        log.error(
            "Science: SUBSTRT1=%d, SUBSTRT2=%d, SUBSIZE1=%d, SUBSIZE2=%d",
            xstart_sci,
            ystart_sci,
            xsize_sci,
            ysize_sci,
        )
        log.error(
            "Reference: SUBSTRT1=%d, SUBSTRT2=%d, SUBSIZE1=%d, SUBSIZE2=%d",
            xstart_ref,
            ystart_ref,
            xsize_ref,
            ysize_ref,
        )
        log.error(
            "Slice indexes: xstart=%d, xstop=%d, ystart=%d, ystop=%d", xstart, xstop, ystart, ystop
        )
        raise ValueError("Bad reference file slice indexes")

    # Extract subarrays from each data attribute in the particular
    # type of reference file model and return a new copy of the
    # data model
    model_type = ref_model.__class__
    sub_model = model_type()
    primary = ref_model.get_primary_array_name()
    for attr in {primary, "err", "dq"}:
        if ref_model.hasattr(attr):
            sub_data = getattr(ref_model, attr)[..., ystart:ystop, xstart:xstop]
            setattr(sub_model, attr, sub_data)
    sub_model.update(ref_model)

    return sub_model


def get_multistripe_subarray_model(sci_model, ref_model):
    """
    Create a multistripe subarray cutout of a reference file.

    Use the multistripe subarray characteristics of a science
    data model to generate a new reference file datamodel,
    containing subarray cutouts of all relevant data arrays
    contained in the reference file model.

    Parameters
    ----------
    sci_model : JWST data model
        The science data model.

    ref_model : JWST data model
        The reference file data model.

    Returns
    -------
    sub_model : JWST data model
        Subarray cutout reference file model.
    """
    # lazy import to avoid circular import error
    from jwst.lib.stripe_utils import stripe_read

    primary = ref_model.get_primary_array_name()
    attrs = {primary}

    for att in ["err", "dq"]:
        if ref_model.hasattr(att):
            attrs.add(att)

    sub_model = stripe_read(sci_model, ref_model, attrs)

    return sub_model


def science_detector_frame_transform(data, fastaxis, slowaxis):
    """
    Swap data array from science to detector frame.

    Use the fastaxis and slowaxis keywords to invert
    and/or transpose data array axes to move from the
    science frame to the detector frame.

    Parameters
    ----------
    data : np.array
        Data array containing at least two dimensions in science orientation.
    fastaxis : int
        Value of the fastaxis keyword for the data array
        to be transformed.
    slowaxis : int
        Value of the slowaxis keyword for the data array
        to be transformed.

    Returns
    -------
    np.array
        Data array transformed to detector orientation.
    """
    # If fastaxis is x-axis
    if np.abs(fastaxis) == 1:
        # Use sign of keywords to possibly reverse the ordering of the axes.
        data = data[..., :: np.sign(slowaxis), :: np.sign(fastaxis)]
    # Else fastaxis is y-axis, also need to transpose array
    else:
        data = data[..., :: np.sign(fastaxis), :: np.sign(slowaxis)]
        data = np.swapaxes(data, -2, -1)
    return data


def detector_science_frame_transform(data, fastaxis, slowaxis):
    """
    Swap data array from detector to science frame.

    Use the fastaxis and slowaxis keywords to invert
    and/or transpose data array axes to move from the
    detector frame to the science frame.

    Parameters
    ----------
    data : np.array
        Data array containing at least two dimensions in detector orientation.
    fastaxis : int
        Value of the fastaxis keyword for the data array
        to be transformed.
    slowaxis : int
        Value of the slowaxis keyword for the data array
        to be transformed.

    Returns
    -------
    np.array
        Data array transformed to science orientation.
    """
    # If fastaxis is x-axis
    if np.abs(fastaxis) == 1:
        # Use sign of keywords to possibly reverse the ordering of the axes.
        data = data[..., :: np.sign(slowaxis), :: np.sign(fastaxis)]
    # Else fastaxis is y-axis, also need to transpose array
    else:
        data = np.swapaxes(data, -2, -1)
        data = data[..., :: np.sign(fastaxis), :: np.sign(slowaxis)]
    return data


class MatchRowError(Exception):
    """Raised when more than one row is matched in a FITS table or list of dict."""

    def __init__(self, message):
        if message is None:
            message = "Expected to match one row only."
        super().__init__(message)


def find_row(ldict, match_keys):
    """
    Find a row in a FITS table matching fields.

    Parameters
    ----------
    ldict : list of dict
        A list of dictionaries, The dictionaries may have any number
        of items but must include all keys in ``match_keys``.
    match_keys : dict
        ``{key: value}`` pairs are matched against all items in ``ldict``
        to find the dict which matches them.

    Returns
    -------
    row : int or None
        FITS table row index, None if no match.

    Raises
    ------
    Warning
        When a field name is not in the table.
    MatchFitsTableRowError
        When more than one rows match.

    Examples
    --------
    >>> ldict = [
    ...     {"row_offset": 2.1, "col_offset": 1.3, "filter": "F444W", "pupil": "CLEAR"},
    ...     {"row_offset": 1, "col_offset": 3, "filter": "F277W", "pupil": "FLAT"},
    ... ]
    >>> match_keys = {"filter": "F444W", "pupil": "CLEAR"}
    >>> result = find_row(ldict, match_keys)
    >>> print(result)
    {'row_offset': 2.1, 'col_offset': 1.3, 'filter': 'F444W', 'pupil': 'CLEAR'}
    """

    def _normalize_strings(field):
        if isinstance(field[0], str):
            return np.array([s.upper() for s in field])
        return field

    # item[1] is always converted to upper case in the `DataSet` initializer.
    results = []
    for d in ldict:
        row = [d[key] == match_keys[key] for key in match_keys]
        if all(row):
            results.append(d)
    if len(results) > 1:
        raise MatchRowError(f"Expected to find one matching row in table, found {len(results)}.")
    if len(results) == 0:
        log.warning("Expected to find one matching row in table, found 0.")
        return None
    return results[0]
