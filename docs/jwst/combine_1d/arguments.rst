Step Arguments
==============
The ``combine_1d`` step has the following optional arguments.

``--exptime_key`` (string, default="exposure_time")
  This is currently not used but kept for future refactoring.

``--sigma_clip`` (float, default=None)
  Optional factor for sigma clipping outliers when combining spectra. If
  a floating point value is provided for ``sigma_clip``, this value will be
  used to set an outlier threshold for any pixels in the input spectra that
  deviate from the median and median absolute deviation of the inputs.
  Defaults to None (such that no clipping is performed).
