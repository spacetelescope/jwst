Step Arguments
==============

The ``pixel_replace`` step has the following step-specific arguments:

``--algorithm`` (str, default="mingrad")
  This sets the method used to estimate flux values for bad pixels. The default is the
  minimum gradient ("mingrad") method, which interpolates values from neighboring
  pixels with the lowest gradient.  Also available is "fit_profile", which uses a profile
  fit to adjacent column values to interpolate missing data. To use only pixels from
  the trace model, set the algorithm to "N/A".

``--use_trace_model`` (bool, default=True)
  Use a trace model attached to the input datamodel to replace pixels, if possible.
  See the :ref:`adaptive_trace_model step <adaptive_trace_model_step>` for more information
  on creating a trace model.

``--n_adjacent_cols`` (int, default=3)
  Number of adjacent columns (on either side of column containing a bad pixel) to use in
  creation of the source profile, in cross-dispersion direction. The total number of
  columns used in the profile will be twice this number; on array edges, the total number
  of columns contributing to the source profile will be less than ``2 * n_adjacent_cols``.
  Ignored when ``algorithm`` is "mingrad" or "N/A".
