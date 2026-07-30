.. _outlier_detection_step_args:

Step Arguments
==============
The ``outlier_detection`` step has the following optional arguments
that control the behavior of the processing.

The following arguments apply to **all modes** unless otherwise specified:

``--snr`` (string, default='5.0 4.0')
  The signal-to-noise values to use for bad pixel identification.
  Since cosmic rays often extend across several pixels, the user
  must specify two cut-off values for determining whether a pixel should
  be masked: the first for detecting the primary cosmic ray, and the
  second (typically lower threshold) for masking lower-level bad pixels
  adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example '5.0 4.0').
  Has no effect for IFU data.

``--save_intermediate_results`` (boolean, default=False)
  Specifies whether or not to save any intermediate products created
  during step processing.

``--good_bits`` (string, default='~DO_NOT_USE')
  The DQ bit values from the input image DQ arrays
  that should be considered 'good'. Any pixel with a DQ value not included
  in this value (or list of values) will be ignored when resampling and flagged
  when building the weight mask. See DQ flag :ref:`dq_parameter_specification` for details.
  Has no effect for IFU data.

The following arguments apply to **imaging or slit-like
spectroscopic data** only:

``--weight_type`` (string, default='ivm')
  The type of data weighting to apply to the resampled data. Available options are 'ivm'
  (default) to compute and use an inverse-variance map, and 'exptime' to
  weight by the exposure time.

``--pixfrac`` (float, default=1.0)
  The pixel fraction used during resampling; valid values go from 0.0 to 1.0.
  Indicates the fraction by which input pixels are "shrunk" before being drizzled onto the
  output image grid. This specifies the size of the footprint, or "dropsize", of a pixel
  in units of the input pixel size.

``--kernel`` (string, default='square')
  The form of the kernel function used to distribute flux onto a
  resampled image. Valid options are 'square' (default), 'point', or 'turbo'.

``--fillval`` (string, default='NAN')
  The value to assign to resampled image pixels that have zero weight or
  do not receive any flux from any input pixels during drizzling.
  Any floating-point value, given as a string, is valid.
  The default value of 'NAN' sets NaN values.

``--maskpt`` (float, default=0.7)
  The percent of maximum weight to use as lower-limit for valid data;
  valid values go from 0.0 to 1.0.

``--scale`` (string, default='1.2 0.7')
  The scaling factor applied to derivative used to identify bad pixels.
  Since cosmic rays often extend across several pixels the user
  must specify two cut-off values for determining whether a pixel should
  be masked: the first for detecting the primary cosmic ray, and the
  second (typically lower threshold) for masking lower-level bad pixels
  adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example '1.2 0.7').

``--backg`` (float, default=0.0)
  User-specified background value to apply to the median image.

``--resample_data`` (boolean, default=True)
  Specifies whether or not to resample the input images when
  performing outlier detection.

``--in_memory`` (boolean, default=True)
  Specifies whether or not to load and create all images that are used during
  processing into memory. If `False`, input files are loaded from disk when
  needed and all intermediate files are stored on disk, rather than in memory.
  Has no effect for spectroscopic data. For imaging data this parameter is
  superseded by the pipeline-level ``in_memory`` parameter set by
  :ref:`calwebb_image3 <calwebb_image3>`.

``--pixmap_stepsize`` (float, default=1.0)
  Indicates the spacing in pixels at which the WCS is evaluated when computing the pixel map.
  Larger step sizes result in faster performance at the cost of accuracy.
  Interpolation is only performed if ``pixmap_stepsize > 1``.
  If it's desired to turn on interpolation, we recommend a value of approx. 10,
  which seemed to work well for most modes during testing.
  Has no effect for spectroscopic data. Has no effect if ``resample_data`` is `False`.
  Default is 1.

``--pixmap_order`` (integer, default=1)
  Interpolating spline order for pixel map computation. Has no effect unless
  ``pixmap_stepsize > 1``. Must be 1 or 3. If it's desired to turn on interpolation,
  we recommend a value of 3, i.e., cubic spline. Default is 1.

The following arguments apply to **IFU data** only:

``--kernel_size`` (string, default='7 7')
  The size of the kernel to use to normalize the pixel differences. The kernel size
  must only contain odd values. Valid values are a pair of integers in a single string
  (for example '7 7', the default).

``--threshold_percent`` (float, default=99.8)
  The threshold (in percent) of the normalized minimum pixel difference used to identify bad pixels.
  Pixels with a normalized minimum difference above this percentage are flagged as outliers.

``--ifu_second_check`` (boolean, default=False)
  Perform a secondary check searching for outliers. This will set outliers
  where ever the difference array of adjacent pixels is a NaN.

The following argument applies to **TSO data** only:

``--rolling_window_width`` (integer, default=25)
  Number of integrations over which to take the median when using rolling-window
  median for TSO observations. The default is 25. If the number of integrations
  is less than or equal to ``rolling_window_width``, a simple median is used instead.
