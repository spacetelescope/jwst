import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def is_subarray(input_model):

    nrows = input_model.data.shape[-2]
    ncols = input_model.data.shape[-1]

    instrument = input_model.meta.instrument.name

    if instrument == 'MIRI':
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


def get_subarray_data(input_model, ref_model):

    # Make sure xstart/ystart exist in science data model
    if (input_model.meta.subarray.xstart is None or
        input_model.meta.subarray.ystart is None):
        raise ValueError('xstart or ystart metadata values ' \
                         + 'not found in input model')

    # Make sure xstart/ystart exist in reference data model
    if (ref_model.meta.subarray.xstart is None or
        ref_model.meta.subarray.ystart is None):
        raise ValueError('xstart or ystart metadata values ' \
                         + 'not found in reference model')

    # Get subarray limits from metadata of input model
    xstart_sci = input_model.meta.subarray.xstart
    xsize_sci = input_model.meta.subarray.xsize
    xstop_sci = xstart_sci + xsize_sci - 1
    ystart_sci = input_model.meta.subarray.ystart
    ysize_sci = input_model.meta.subarray.ysize
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

    # Compute slice limits, in 1-indexed units
    xstart = xstart_sci - xstart_ref + 1
    xstop = xstart + xsize_sci - 1
    ystart = ystart_sci - ystart_ref + 1
    ystop = ystart + ysize_sci - 1
    log.debug('slice xstart=%d, xstop=%d, ystart=%d, ystop=%d',
              xstart, xstop, ystart, ystop)

    if (xstart < 1 or ystart < 1 or
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

    return ref_model.data[ystart-1:ystop, xstart-1:xstop]
