Step Arguments
==============

The ``pixel_replace`` step has the following step-specific arguments:

``--algorithm`` (str, default='fit_profile')
  This sets the method used to estimate flux values for bad pixels. The default 'fit_profile' uses a profile
  fit to adjacent column values.  The minimum gradient ('mingrad') method is also available for the MIRI MRS.

``--n_adjacent_cols`` (int, default=3)
  Number of adjacent columns (on either side of column containing a bad pixel) to use in
  creation of the source profile, in cross-dispersion direction. The total number of
  columns used in the profile will be twice this number; on array edges, the total number
  of columns contributing to the source profile will be less than ``2 * n_adjacent_cols``.

