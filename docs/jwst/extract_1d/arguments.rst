Step Arguments
==============

The extract_1d step has five step-specific arguments.

*  ``--smoothing_length``

If ``smoothing_length`` is greater than 1 (and is an odd integer), the
background will be smoothed in the dispersion direction with a boxcar of
this width.  If ``smoothing_length`` is None (the default), the step will
attempt to read the value from the reference file.  If a value was
specified in the reference file, that will be used.  Note that in this
case a different value can be specified for each slit.  If no value was
specified either by the user or in the reference file, no background
smoothing will be done.

*  ``--bkg_order``

This is the order of a polynomial function to be fit to the background
regions.  The fit is done independently for each column (or row, if the
dispersion is vertical) of the input image, and the fitted curve will be
subtracted from the target data.  ``bkg_order`` = 0 (the minimum allowed
value) means to fit a constant.  The user-supplied value (if any)
overrides the value in the reference file.  If neither is specified, a
value of 0 will be used.

*  ``--log_increment``

Most log messages are suppressed while looping over integrations, i.e. when
the input is a CubeModel or a 3-D SlitModel.  Messages will be logged while
processing the first integration, but since they would be the same for
every integration, most messages will only be written once.  However, since
there can be hundreds or thousands of integrations, which can take a long
time to process, it would be useful to log a message every now and then to
let the user know that the step is still running.

``log_increment`` is an integer, with default value 50.  If it is greater
than 0, an INFO message will be printed every ``log_increment``
integrations, e.g. "... 150 integrations done".

*  ``--subtract_background``

This is a boolean flag to specify whether the background should be
subtracted.  If None, the value in the extract_1d reference file (if any)
will be used.  If not None, this parameter overrides the value in the
extract_1d reference file.

*  ``--apply_nod_offset``

This is a boolean flag to specify whether the target and background positions
specified in the reference file should be shifted to account for nod
and/or dither offset.  If None (the default), the value in the reference
file will be used, or it will be set to True if it is not specified in
the reference file.  The offset is determined by finding the location in
the data corresponding to the target position (as given by keywords
TARG_RA and TARG_DEC).

At the time of writing, a nod/dither offset will not be applied if the
source is extended.  It will also not be applied for wide-field slitless
spectroscopy data, or NIRSpec fixed-slit, or NIRSpec MOS (MSA) data.
