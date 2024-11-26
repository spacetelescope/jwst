.. _outlier_detection_step_args:

Step Arguments
==============

The outlier detection step has the following optional arguments
that control the behavior of the processing.
For more details about step arguments (including datatypes, possible values
and defaults) see :py:obj:`jwst.outlier_detection.OutlierDetectionStep.spec`.


General Step Arguments
----------------------
The following arguments apply to all modes unless otherwise specified:

``--save_intermediate_results``
  Specifies whether or not to save any intermediate products created
  during step processing.

``--good_bits``
  The DQ bit values from the input image DQ arrays
  that should be considered 'good'. Any pixel with a DQ value not included
  in this value (or list of values) will be ignored when resampling and flagged
  when building the weight mask. See DQ flag :ref:`dq_parameter_specification` for details.
  Has no effect for IFU data.

``--snr``
  The signal-to-noise values to use for bad pixel identification.
  Since cosmic rays often extend across several pixels, the user
  must specify two cut-off values for determining whether a pixel should
  be masked: the first for detecting the primary cosmic ray, and the
  second (typically lower threshold) for masking lower-level bad pixels
  adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example "5.0 4.0").
  Has no effect for IFU data.


Step Arguments for Imaging and Slit-like Spectroscopic data
-----------------------------------------------------------

``--weight_type``
  The type of data weighting to apply to the resampled data. Available options are ``ivm``
  (default) to compute and use an inverse-variance map, and ``exptime`` to
  weight by the exposure time.

``--pixfrac``
  The pixel fraction used during resampling; valid values go from 0.0 to 1.0.
  Indicates the fraction by which input pixels are "shrunk" before being drizzled onto the
  output image grid. This specifies the size of the footprint, or "dropsize", of a pixel
  in units of the input pixel size.

``--kernel``
  The form of the kernel function used to distribute flux onto a
  resampled image.

``--fillval``
  The value to assign to resampled image pixels that have zero weight or
  do not receive any flux from any input pixels during drizzling.
  Any floating-point value, given as a string, is valid.
  The default value of 'NAN' sets NaN values.

``--maskpt``
  The percent of maximum weight to use as lower-limit for valid data;
  valid values go from 0.0 to 1.0.

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

``--resample_data``
  Specifies whether or not to resample the input images when
  performing outlier detection.

``--in_memory``
  Specifies whether or not to load and create all images that are used during
  processing into memory. If ``False``, input files are loaded from disk when
  needed and all intermediate files are stored on disk, rather than in memory.
  Has no effect for spectroscopic data. For imaging data this parameter is 
  superseded by the pipeline-level ``in_memory`` parameter set by
  ``calwebb_image3``.


Step Arguments for IFU data
---------------------------

``--kernel_size``
  The size of the kernel to use to normalize the pixel differences. The kernel size
  must only contain odd values. Valid values are a pair of ints in a single string
  (for example "7 7", the default).

``--threshold_percent``
  The threshold (in percent) of the normalized minimum pixel difference used to identify bad pixels.
  Pixels with a normalized minimum difference above this percentage are flagged as outliers.

``--ifu_second_check``
  Perform a secondary check searching for outliers. This will set outliers
  where ever the difference array of adjacent pixels is a Nan.


Step Arguments for TSO data
---------------------------

``--rolling_window_width``
  Number of integrations over which to take the median when using rolling-window
  median for TSO observations. The default is 25. If the number of integrations
  is less than or equal to ``rolling_window_width``, a simple median is used instead.


Step Arguments for Coronagraphic data
-------------------------------------
General step arguments apply to coronagraphic data. No additional arguments are used.