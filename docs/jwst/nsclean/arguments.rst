Step Arguments
==============

The ``nsclean`` step has the following optional arguments to control
the behavior of the processing.

``--n_sigma`` (float, default=5.0)
  The sigma-clipping threshold to use when searching for outliers
  and illuminated pixels to be excluded from use in the fitting
  process.

``--save_mask`` (boolean, default=False)
  A flag to indicate whether the mask constructed by the step
  should be saved to a file.
