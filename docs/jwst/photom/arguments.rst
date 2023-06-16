Arguments
=========
The ``photom`` step has the following optional arguments.

``--inverse`` (boolean, default=False)
  A flag to indicate whether the math operations used to apply the
  correction should be inverted (i.e. divide the calibration data
  into the science data, instead of the usual multiplication).

``--source_type`` (string, default=None)
  Force the processing to use the given source type (POINT, EXTENDED),
  instead of using the information contained in the input data.

``--mrs_time_correction`` (boolean, default=True)
   A flag to indicate whether to turn on the time and wavelength dependent
   correction for MIRI MRS data. 
   
