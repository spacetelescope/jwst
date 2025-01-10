.. _outlier_detection_coron_step_args:

Step Arguments
==============

The outlier detection step for coronagraphic data has the following
optional arguments that control the behavior of the processing.
For more details about step arguments (including datatypes, possible values
and defaults) see :py:obj:`jwst.outlier_detection.OutlierDetectionCoronStep.spec`.

``--save_intermediate_results``
  Specifies whether or not to save any intermediate products created
  during step processing.

``--good_bits``
  The DQ bit values from the input image DQ arrays
  that should be considered 'good'. Any pixel with a DQ value not included
  in this value (or list of values) will be ignored when resampling and flagged
  when building the weight mask. See DQ flag :ref:`dq_parameter_specification` for details.

``--snr``
  The signal-to-noise values to use for bad pixel identification.
  Since cosmic rays often extend across several pixels, the user
  must specify two cut-off values for determining whether a pixel should
  be masked: the first for detecting the primary cosmic ray, and the
  second (typically lower threshold) for masking lower-level bad pixels
  adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example "5.0 4.0").

``--maskpt``
  The percent of maximum weight to use as lower-limit for valid data;
  valid values go from 0.0 to 1.0.
