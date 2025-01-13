.. _outlier_detection_ifu_step_args:

Step Arguments
==============

The outlier detection step for IFU data has the following
optional arguments that control the behavior of the processing.
For more details about step arguments (including datatypes, possible values
and defaults) see :py:obj:`jwst.outlier_detection.OutlierDetectionIFUStep.spec`.

``--save_intermediate_results``
  Specifies whether or not to save any intermediate products created
  during step processing.

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
