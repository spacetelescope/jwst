Step Arguments
==============

The ``pixel_replace`` step has the following step-specific arguments:

``--algorithm`` (str, default="mingrad")
  This sets the method used to estimate flux values for bad pixels. The default is the
  minimum gradient ("mingrad") method, which interpolates values from neighboring
  pixels with the lowest gradient. Also available is "fit_profile", which uses a simple profile
  fit to adjacent column values to interpolate missing data, and "trace_model" which
  uses a detailed spectral trace model to replace pixels, if possible.

  See the :ref:`adaptive_trace_model step <adaptive_trace_model_step>` for more information
  on creating a trace model. If that step has not been run prior
  to calling ``pixel_replace`` with ``algorithm=trace_model``, then it will be run on all input
  data at the beginning of the step, with ``oversample=1`` and otherwise default parameters.

``--n_adjacent_cols`` (int, default=3)
  Number of adjacent columns (on either side of column containing a bad pixel) to use in
  creation of the source profile, in cross-dispersion direction. The total number of
  columns used in the profile will be twice this number; on array edges, the total number
  of columns contributing to the source profile will be less than ``2 * n_adjacent_cols``.
  Ignored when ``algorithm`` is "mingrad" or "trace_model".
