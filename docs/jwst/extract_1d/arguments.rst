Step Arguments
==============

The extract_1d step has two step-specific arguments.  Currently neither of
these is used for IFU data.

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
