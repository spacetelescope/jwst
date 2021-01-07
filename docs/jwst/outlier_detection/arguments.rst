.. _outlier_detection_step_args:

Step Arguments
==============
The `outlier_detection` step has the following optional arguments
that control the behavior of the processing:

``--weight_type`` (string, default='exptime')
  The type of data weighting to use during resampling;
  options are 'exptime', 'error', and 'None'.

``--pixfrac`` (float, default=1.0)
  The pixel fraction used during resampling;
  valid values go from 0.0 to 1.0.

``--kernel`` (string, default='square')
  The form of the kernel function used to distribute flux onto a
  resampled image. Options are 'square', 'turbo', 'point',
  'lanczos', and 'tophat'.

``--fillval`` (string, default='INDEF')
  The value to assign to resampled image pixels that have zero weight or
  do not receive any flux from any input pixels during drizzling.
  Any floating-point value, given as a string, is valid.
  A value of 'INDEF' will use the last zero weight flux.

``--nlow`` (integer, default=0)
  The number of low values in each pixel stack to ignore
  when computing the median value.

``--nhigh`` (integer, default=0)
  The number of high values in each pixel stack to ignore
  when computing the median value.

``--maskpt`` (float, default=0.7)
  The percent of maximum weight to use as lower-limit for valid data;
  valid values go from 0.0 to 1.0.

``--grow`` (integer, default=1)
  The radius, in pixels, from a bad pixel for neighbor rejection.

``--snr`` (string, default='4.0 3.0')
  The signal-to-noise values to use for bad pixel identification. Valid
  values are a pair of floating-point values in a single string.

``--scale`` (string, default='0.5 0.4')
  The scaling factor applied to derivative used to identify bad pixels.
  Valid values are a pair of floating-point values in a single string.

``--backg`` (float, default=0.0)
  User-specified background value to apply to the median image.

``--save_intermediate_results`` (boolean, default=False)
  Specifies whether or not to save any intermediate products created
  during step processing.

``--resample_data`` (boolean, default=True)
  Specifies whether or not to resample the input images when
  performing outlier detection.

``--good_bits`` (string, default="~DO_NOT_USE")
  The DQ bit values from the input image DQ arrays
  that should be considered 'good' when building the weight mask. See
  DQ flag :ref:`dq_parameter_specification` for details.

``--scale_detection`` (bool, default=False)
  Specifies whether or not to rescale the individual input images
  to match total signal when doing comparisons.

``--allowed_memory`` (float, default=None)
  Specifies the fractional amount of
  free memory to allow when creating the resampled image. If ``None``, the
  environmental variable ``DMODEL_ALLOWED_MEMORY`` is used. If not defined, no
  check is made. If the resampled image would be larger than specified, an
  ``OutputTooLargeError`` exception will be generated.

  For example, if set to ``0.5``, only resampled images that use less than half
  the available memory can be created.
