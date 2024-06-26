.. _outlier_detection_step_args:

For more details about step arguments (including datatypes, possible values
and defaults) see :py:obj:`jwst.outlier_detection.OutlierDetectionStep.spec`.

Step Arguments for Non-IFU data
===============================
The `outlier_detection` step for non-IFU data has the following optional arguments
that control the behavior of the processing:

``--weight_type``
  The type of data weighting to use during resampling.

``--pixfrac``
  The pixel fraction used during resampling;
  valid values go from 0.0 to 1.0.

``--kernel``
  The form of the kernel function used to distribute flux onto a
  resampled image.

``--fillval``
  The value to assign to resampled image pixels that have zero weight or
  do not receive any flux from any input pixels during drizzling.
  Any floating-point value, given as a string, is valid.
  A value of 'INDEF' will use the last zero weight flux.

``--nlow``
  Deprecated and has no effect. This parameter will be removed
  in a future version.

``--nhigh``
  Deprecated and has no effect. This parameter will be removed
  in a future version.

``--maskpt``
  The percent of maximum weight to use as lower-limit for valid data;
  valid values go from 0.0 to 1.0.

``--snr``
  The signal-to-noise values to use for bad pixel identification.
  Since cosmic rays often extend across several pixels the user
  must specify two cut-off values for determining whether a pixel should
  be masked: the first for detecting the primary cosmic ray, and the
  second (typically lower threshold) for masking lower-level bad pixels
  adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example "5.0 4.0").

``--scale``
  The scaling factor applied to derivative used to identify bad pixels.
  Since cosmic rays often extend across several pixels the user
  must specify two cut-off values for determining whether a pixel should
  be masked: the first for detecting the primary cosmic ray, and the
  second (typically lower threshold) for masking lower-level bad pixels
  adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example "1.2 0.7").

``--backg``
  User-specified background value to apply to the median image.

``--rolling_window_width``
  Number of integrations over which to take the median when using rolling-window
  median for TSO observations.

``--save_intermediate_results``
  Specifies whether or not to save any intermediate products created
  during step processing.

``--resample_data``
  Specifies whether or not to resample the input images when
  performing outlier detection.

``--good_bits``
  The DQ bit values from the input image DQ arrays
  that should be considered 'good' when building the weight mask. See
  DQ flag :ref:`dq_parameter_specification` for details.

``--allowed_memory``
  Specifies the fractional amount of
  free memory to allow when creating the resampled image. If ``None``, the
  environment variable ``DMODEL_ALLOWED_MEMORY`` is used. If not defined, no
  check is made. If the resampled image would be larger than specified, an
  ``OutputTooLargeError`` exception will be generated.

  For example, if set to ``0.5``, only resampled images that use less than half
  the available memory can be created.

``--in_memory``
  Specifies whether or not to load and create all images that are used during
  processing into memory. If ``False``, input files are loaded from disk when
  needed and all intermediate files are stored on disk, rather than in memory.

Step Arguments for IFU data
===========================
The `outlier_detection` step for IFU data has the following optional arguments
that control the behavior of the processing:

``--kernel_size``
  The size of the kernel to use to normalize the pixel differences. The kernel size
  must only contain odd values. Valid values are a pair of ints in a single string
  (for example "7 7").

``--threshold_percent``
  The threshold (in percent) of the normalized minimum pixel difference used to identify bad pixels.
  Pixels with   a normalized minimum pixel difference above this percentage are flagged as a outlier.

``--save_intermediate_results``
  Specifies whether or not to save any intermediate products created
  during step processing.

``--in_memory``
  Specifies whether or not to load and create all images that are used during
  processing into memory. If ``False``, input files are loaded from disk when
  needed and all intermediate files are stored on disk, rather than in memory.
