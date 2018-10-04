from jwst import datamodels
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def is_subarray(input_model):

    """
    Check to see if a data model comes from a subarray readout.
    If the data dimensions are less than 2048x2048 (or 1032x1024
    for MIRI), it is assumed to be a subarray.

    Parameters
    ----------
    input_model: JWST data model
        input data model to be checked

    Returns
    -------
    True is the primary data array of the model is smaller than
    full-frame, False otherwise.
    """

    # Get the dimensions of the model's primary data array
    nrows = input_model.data.shape[-2]
    ncols = input_model.data.shape[-1]

    # Compare the dimensions to a full-frame
    if input_model.meta.instrument.name.upper() == 'MIRI':
        if ncols < 1032 or nrows < 1024:
            log.debug('Input exposure is a subarray readout')
            return True
        else:
            return False
    else:
        if ncols < 2048 or nrows < 2048:
            log.debug('Input exposure is a subarray readout')
            return True
        else:
            return False


def ref_matches_sci(sci_model, ref_model):

    """
    Short Summary
    -------------
    Check to see if the science model has the same subarray
    characteristics as the reference model. Also performs error
    checking on the subarray keywords in both the reference file
    and science data.

    Parameters
    ----------
    sci_model: JWST data model
        science data model

    ref_model: JWST data model
        reference file data model

    Returns
    -------
    True if the science model has the same subarray parameters
    as the reference model, False otherwise.
    """

    # Get the reference file subarray parameters
    xstart_ref = ref_model.meta.subarray.xstart
    xsize_ref = ref_model.meta.subarray.xsize
    ystart_ref = ref_model.meta.subarray.ystart
    ysize_ref = ref_model.meta.subarray.ysize

    # Make sure the attributes were populated
    if (xstart_ref is None) or (xsize_ref is None) or \
       (ystart_ref is None) or (ysize_ref is None):

        # If the ref file is full-frame, set the missing params to
        # default values
        if ref_model.meta.instrument.name.upper() == 'MIRI':
            if (ref_model.data.shape[-1] == 1032 and \
                ref_model.data.shape[-2] == 1024):
                log.warning('Missing subarray corner/size keywords in reference file')
                log.warning('Setting them to full-frame default values')
                xstart_ref = 1
                ystart_ref = 1
                xsize_ref = 1032
                ysize_ref = 1024
                ref_model.meta.subarray.xstart = xstart_ref
                ref_model.meta.subarray.ystart = ystart_ref
                ref_model.meta.subarray.xsize = xsize_ref
                ref_model.meta.subarray.ysize = ysize_ref
            else:
                log.error('Missing subarray corner/size keywords in reference file')
                raise ValueError("Can't determine ref file subarray properties")
        else:
            if (ref_model.data.shape[-1] == 2048) and \
               (ref_model.data.shape[-2] == 2048):
                log.warning('Missing subarray corner/size keywords in reference file')
                log.warning('Setting them to full-frame default values')
                xstart_ref = 1
                ystart_ref = 1
                xsize_ref = 2048
                ysize_ref = 2048
                ref_model.meta.subarray.xstart = xstart_ref
                ref_model.meta.subarray.ystart = ystart_ref
                ref_model.meta.subarray.xsize = xsize_ref
                ref_model.meta.subarray.ysize = ysize_ref
            else:
                log.error('Missing subarray corner/size keywords in reference file')
                raise ValueError("Can't determine ref file subarray properties")

    log.debug(" ref xstart=%d, xsize=%d, ystart=%d, ysize=%d" % \
        (xstart_ref, xsize_ref, ystart_ref, ysize_ref))

    # Make sure the starting corners are valid
    if (xstart_ref < 1) or (ystart_ref < 1):
        log.error("Reference file subarray corners are invalid")
        raise ValueError("Bad subarray corners")

    # Make sure the size attributes are consistent with data array
    if hasattr(ref_model, 'coeffs'):
        xsize_data = ref_model.coeffs.shape[-1]
        ysize_data = ref_model.coeffs.shape[-2]
    else:
        xsize_data = ref_model.shape[-1]
        ysize_data = ref_model.shape[-2]

    if (xsize_ref != xsize_data) or (ysize_ref != ysize_data):
        # Check for NIRSpec IRS2 mode, which is allowed to have a mismatch
        if (ysize_ref == 2048) and (ysize_data == 3200):
            ysize_ref = ysize_data
        else:
            log.warning("Reference file data array size doesn't match subarray params")
            log.warning("Using actual array size")
            xsize_ref = xsize_data
            ysize_ref = ysize_data

    # Get the science model subarray parameters
    try:
        # If the science data model has already gone through
        # extract_2d, the metadata is in the top level of the model
        xstart_sci = sci_model.xstart
        xsize_sci = sci_model.xsize
        ystart_sci = sci_model.ystart
        ysize_sci = sci_model.ysize
    except:
        # Otherwise the metadata is in the meta tree
        xstart_sci = sci_model.meta.subarray.xstart
        xsize_sci = sci_model.meta.subarray.xsize
        ystart_sci = sci_model.meta.subarray.ystart
        ysize_sci = sci_model.meta.subarray.ysize

    # Make sure the attributes were populated
    if (xstart_sci is None) or (xsize_sci is None) or \
       (ystart_sci is None) or (ysize_sci is None):
        log.error('Missing subarray corner/size keywords in science file')
        raise ValueError("Can't determine science file subarray properties")
    else:
        log.debug(" sci xstart=%d, xsize=%d, ystart=%d, ysize=%d" % \
            (xstart_sci, xsize_sci, ystart_sci, ysize_sci))

    # Make sure the starting corners are valid
    if (xstart_sci < 1) or (ystart_sci < 1):
        log.error("Science file subarray corners are invalid")
        raise ValueError("Bad subarray corners")

    # Make sure the size attributes are consistent with data array
    if (xsize_sci != sci_model.shape[-1]) or \
       (ysize_sci != sci_model.shape[-2]):

        # Check for NIRSpec IRS2 mode, which is allowed to have a mismatch
        if (ysize_sci == 2048) and (sci_model.shape[-2] == 3200):
            ysize_sci = sci_model.shape[-2]
        else:
            log.warning("Science file data array size doesn't match subarray params")
            log.warning("Using actual array size")
            xsize_sci = sci_model.shape[-1]
            ysize_sci = sci_model.shape[-2]

    # Finally, see if all of the science file subarray params match those
    # of the reference file
    if (xstart_ref == xstart_sci and xsize_ref == xsize_sci and
        ystart_ref == ystart_sci and ysize_ref == ysize_sci):
        return True
    else:
        return False


def get_subarray_data(sci_model, ref_model):

    """
    Extract a subarray from the data attribute of a reference file
    data model that matches the subarray characteristics of a
    science data model. Only the extracted data array is returned.

    Parameters
    ----------
    sci_model: JWST data model
        science data model

    ref_model: JWST data model
        reference file data model

    Returns
    -------
    array: 2-D extracted data array
    """

    # Make sure xstart/ystart exist in science data model
    if (sci_model.meta.subarray.xstart is None or
        sci_model.meta.subarray.ystart is None):

        # If the science file is full-frame, set the missing params to
        # default values
        if (sci_model.data.shape[-1] == 2048) and \
           (sci_model.data.shape[-2] == 2048):
            sci_model.meta.subarray.xstart = 1
            sci_model.meta.subarray.ystart = 1
            sci_model.meta.subarray.xsize = 2048
            sci_model.meta.subarray.ysize = 2048
        else:
            raise ValueError('xstart or ystart metadata values ' \
                             + 'not found in input model')

    # Make sure xstart/ystart exist in reference data model
    if (ref_model.meta.subarray.xstart is None or
        ref_model.meta.subarray.ystart is None):

        # If the ref file is full-frame, set the missing params to
        # default values
        if ref_model.meta.instrument.name.upper() == 'MIRI':
            if (ref_model.data.shape[-1] == 1032 and \
                ref_model.data.shape[-2] == 1024):
                ref_model.meta.subarray.xstart = 1
                ref_model.meta.subarray.ystart = 1
                ref_model.meta.subarray.xsize = 1032
                ref_model.meta.subarray.ysize = 1024
            else:
                raise ValueError('xstart or ystart metadata values ' \
                                 + 'not found in reference model')
        else:
            if (ref_model.data.shape[-1] == 2048) and \
               (ref_model.data.shape[-2] == 2048):
                ref_model.meta.subarray.xstart = 1
                ref_model.meta.subarray.ystart = 1
                ref_model.meta.subarray.xsize = 2048
                ref_model.meta.subarray.ysize = 2048
            else:
                raise ValueError('xstart or ystart metadata values ' \
                                 + 'not found in reference model')

    # Get subarray limits from metadata of input model
    xstart_sci = sci_model.meta.subarray.xstart
    xsize_sci = sci_model.meta.subarray.xsize
    xstop_sci = xstart_sci + xsize_sci - 1
    ystart_sci = sci_model.meta.subarray.ystart
    ysize_sci = sci_model.meta.subarray.ysize
    ystop_sci = ystart_sci + ysize_sci - 1
    log.debug('science xstart=%d, xstop=%d, ystart=%d, ystop=%d',
              xstart_sci, xstop_sci, ystart_sci, ystop_sci)

    # Get subarray limits from metadata of reference model
    xstart_ref = ref_model.meta.subarray.xstart
    xsize_ref = ref_model.meta.subarray.xsize
    xstop_ref = xstart_ref + xsize_ref - 1
    ystart_ref = ref_model.meta.subarray.ystart
    ysize_ref = ref_model.meta.subarray.ysize
    ystop_ref = ystart_ref + ysize_ref - 1
    log.debug('reference xstart=%d, xstop=%d, ystart=%d, ystop=%d',
              xstart_ref, xstop_ref, ystart_ref, ystop_ref)

    # Compute slice limits, in 0-indexed python notation
    xstart = xstart_sci - xstart_ref
    ystart = ystart_sci - ystart_ref
    xstop = xstart + xsize_sci
    ystop = ystart + ysize_sci
    log.debug('slice xstart=%d, xstop=%d, ystart=%d, ystop=%d',
              xstart, xstop, ystart, ystop)

    # Make sure that the slice limits are within the bounds of
    # the reference file data array
    if (xstart < 0 or ystart < 0 or
        xstop > ref_model.meta.subarray.xsize or
        ystop > ref_model.meta.subarray.ysize):
        log.error('Computed reference file slice indexes are ' \
                  + 'incompatible with size of reference data array')
        log.error('xstart=%d, xstop=%d, ystart=%d, ystop=%d',
                   xstart, xstop, ystart, ystop)
        log.error('Reference xsize=%d, ysize=%d',
                   ref_model.meta.subarray.xsize,
                   ref_model.meta.subarray.ysize)
        raise ValueError('Bad reference file slice indexes')

    return ref_model.data[ystart:ystop, xstart:xstop]


def get_subarray_model(sci_model, ref_model):

    """
    Create a subarray version of a reference file model that matches
    the subarray characteristics of a science data model. An new
    model is created that contains subarrays of all data arrays
    contained in the reference file model.

    Parameters
    ----------
    sci_model: JWST data model
        science data model

    ref_model: JWST data model
        reference file data model

    Returns
    -------
    sub_model: JWST data model
        subarray version of the reference file model
    """

    # Get the science model subarray params
    xstart_sci = sci_model.meta.subarray.xstart
    xsize_sci = sci_model.meta.subarray.xsize
    ystart_sci = sci_model.meta.subarray.ystart
    ysize_sci = sci_model.meta.subarray.ysize

    # Get the reference model subarray params
    xstart_ref = ref_model.meta.subarray.xstart
    xsize_ref = ref_model.meta.subarray.xsize
    ystart_ref = ref_model.meta.subarray.ystart
    ysize_ref = ref_model.meta.subarray.ysize

    # Compute the slice indexes, in 0-indexed python frame
    xstart = xstart_sci - xstart_ref
    ystart = ystart_sci - ystart_ref
    xstop = xstart + xsize_sci
    ystop = ystart + ysize_sci
    log.debug("slice xstart=%d, xstop=%d, ystart=%d, ystop=%d",
              xstart, xstop, ystart, ystop)

    # Make sure that the slice limits are within the bounds of
    # the reference file data array
    if (xstart < 0 or ystart < 0 or
        xstop > ref_model.meta.subarray.xsize or
        ystop > ref_model.meta.subarray.ysize):
        log.error('Computed reference file slice indexes are ' \
                  + 'incompatible with size of reference data array')
        log.error('xstart=%d, xstop=%d, ystart=%d, ystop=%d',
                   xstart, xstop, ystart, ystop)
        log.error('Reference xsize=%d, ysize=%d',
                   ref_model.meta.subarray.xsize,
                   ref_model.meta.subarray.ysize)
        raise ValueError('Bad reference file slice indexes')

    # Extract subarrays from each data attribute in the particular
    # type of reference file model and return a new copy of the
    # data model
    if isinstance(ref_model, datamodels.FlatModel):
        sub_data = ref_model.data[ystart:ystop, xstart:xstop]
        sub_err = ref_model.err[ystart:ystop, xstart:xstop]
        sub_dq = ref_model.dq[ystart:ystop, xstart:xstop]
        sub_model = datamodels.FlatModel(data=sub_data, err=sub_err, dq=sub_dq)
        sub_model.update(ref_model)
    elif isinstance(ref_model, datamodels.GainModel):
        sub_data = ref_model.data[ystart:ystop, xstart:xstop]
        sub_model = datamodels.GainModel(data=sub_data)
        sub_model.update(ref_model)
    elif isinstance(ref_model, datamodels.LinearityModel):
        sub_data = ref_model.coeffs[:, ystart:ystop, xstart:xstop]
        sub_dq = ref_model.dq[ystart:ystop, xstart:xstop]
        sub_model = datamodels.LinearityModel(coeffs=sub_data, dq=sub_dq)
        sub_model.update(ref_model)
    elif isinstance(ref_model, datamodels.MaskModel):
        sub_dq = ref_model.dq[ystart:ystop, xstart:xstop]
        sub_model = datamodels.MaskModel(dq=sub_dq)
        sub_model.update(ref_model)
    elif isinstance(ref_model, datamodels.ReadnoiseModel):
        sub_data = ref_model.data[ystart:ystop, xstart:xstop]
        sub_model = datamodels.ReadnoiseModel(data=sub_data)
        sub_model.update(ref_model)
    elif isinstance(ref_model, datamodels.SaturationModel):
        sub_data = ref_model.data[ystart:ystop, xstart:xstop]
        sub_dq = ref_model.dq[ystart:ystop, xstart:xstop]
        sub_model = datamodels.SaturationModel(data=sub_data, dq=sub_dq)
        sub_model.update(ref_model)
    elif isinstance(ref_model, datamodels.SuperBiasModel):
        sub_data = ref_model.data[ystart:ystop, xstart:xstop]
        sub_err = ref_model.err[ystart:ystop, xstart:xstop]
        sub_dq = ref_model.dq[ystart:ystop, xstart:xstop]
        sub_model = datamodels.SuperBiasModel(data=sub_data, err=sub_err, dq=sub_dq)
        sub_model.update(ref_model)
    else:
        log.warning('Unsupported reference file model type')
        sub_model = None

    return sub_model
