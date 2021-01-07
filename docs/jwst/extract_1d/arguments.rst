Step Arguments
==============

The ``extract_1d`` step has the following step-specific arguments.

``--smoothing_length``
  If ``smoothing_length`` is greater than 1 (and is an odd integer), the
  image data used to perform background extraction will be smoothed in the
  dispersion direction with a boxcar of this width.  If ``smoothing_length``
  is None (the default), the step will attempt to read the value from the
  EXTRACT1D reference file.  If a value is specified in the reference file,
  that value will be used.  Note that by specifying this parameter in the
  EXTRACT1D reference file a different value can be designated for each slit.
  If no value is specified either by the user or in the EXTRACT1D reference
  file, no background smoothing is done.

``--bkg_fit``
  The type of fit to perform to the background data in each image column
  (or row, if the dispersion is vertical). There are three allowed values:
  "poly" (the default), "mean", and "median". If set to "poly", the background
  values for each pixel within all background regions in a given column (or
  row) will be fit with a polynomial of order "bkg_order" (see below).
  Values of "mean" and "median" compute the simple average and median,
  respectively, of the background region values in each column (or row).
  This parameter can also be specified in the EXTRACT1D reference file. If
  specified in the reference file and given as an argument to the step by
  the user, the user-supplied value takes precedence.

``--bkg_order``
  The order of a polynomial function to be fit to the background
  regions.  The fit is done independently for each column (or row, if the
  dispersion is vertical) of the input image, and the fitted curve will be
  subtracted from the target data.  ``bkg_order`` = 0 (the minimum allowed
  value) means to fit a constant.  The user-supplied value (if any)
  overrides the value in the EXTRACT1D reference file.  If neither is specified, a
  value of 0 will be used. If a sufficient number of valid data points is
  unavailable to construct the polynomial fit, the fit will be forced to
  0 for that particular column (or row). If "bkg_fit" is not "poly", this
  parameter will be ignored.

``--log_increment``
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

``--subtract_background``
  This is a boolean flag to specify whether the background should be
  subtracted.  If None, the value in the EXTRACT1D reference file (if any)
  will be used.  If not None, this parameter overrides the value in the
  reference file.

``--use_source_posn``
  This is a boolean flag to specify whether the target and background extraction
  region locations specified in the EXTRACT1D reference file should be shifted
  to account for the expected position of the source. If None (the default),
  the step will make the decision of whether to use the source position based
  on the observing mode and the source type. The source position will only be
  used for point sources and for modes where the source could be located
  off-center due to things like nodding or dithering. If turned on, the sky
  (RA/Dec) position of the source is used in conjunction with the World
  Coordinate System (WCS) to compute the x/y source location. For long-slit
  type modes (e.g. MIRI LRS and NIRSpec fixed-slit and MOS), only the position
  in the cross-dispersion direction is used to potentially offset the
  extraction regions in that direction.

``--apply_apcorr``
  Switch to select whether or not to apply an APERTURE correction during the
  Extract1dStep processing. Default is ``True``
