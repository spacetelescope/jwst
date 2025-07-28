.. _msb_step_args:

Step Arguments
==============
The master background subtraction step uses the following optional arguments.

``--median_kernel``
  Optional user-supplied kernel size for a moving-median boxcar filter to filter 
  outliers in the master background spectrum.  The kernel size must be an odd integer.
  Even integers will be rounded down to the nearest odd integer.
  Defaults to 1 (which applies no median filtering).

``--user_background``
  The file name of a user-supplied 1-D master background spectrum. Must be in the form
  of a standard :ref:`x1d <x1d>` product containing a single 'EXTRACT1D' extension.
  When a user background spectrum is supplied, it is used for the subtraction instead of
  a computed master background, and the name of the user-supplied file is recorded in the
  MSTRBKGD keyword in the output product(s).
  Defaults to ``None``.

``--save_background``
  A boolean indicating whether the computed 1-D master background spectrum should be saved
  to a file. The file name uses a product type suffix of "masterbg".
  For the `master_background_mos` step, multiple files will be produced including the 1-D 
  master background spectrum (saved with the suffix "masterbg1d"), the expanded 2-D background spectra
  for each MOS slitlet (with the suffix "masterbg2d"), and the 1-D background spectra 
  that were combined into the master background spectrum (with the suffix "bkgx1d").
  If a user-supplied background is specified, this argument is ignored.
  Defaults to ``False``.

``--force_subtract``
  A boolean indicating whether or not to override the step's built-in logic for determining
  if the step should be applied. By default, the step will be skipped if the
  :ref:`calwebb_spec2 <calwebb_spec2>` :ref:`background <background_subtraction>` step has
  already been applied. If ``--force_subtract = True``, the master background will be
  applied.

``--output_use_model``
  A boolean indicating whether to use the "filename" meta attribute in the data model to
  determine the name of the output file created by the step. Defaults to ``True``.

``--sigma_clip``
  Factor for sigma clipping outliers and contaminated spectra when combining MOS 
  background spectra in the `master_background_mos` step.  The value of ``sigma_clip`` 
  will be used to set an outlier threshold for clipping any pixels in the background 
  spectra that deviate from the median and median absolute deviation of the inputs before
  combining the background spectra.  Setting ``sigma_clip`` to None will
  skip any outlier clipping.  This parameter is only available in `master_background_mos`
  step and is not available in the generic `master_background` step.
  Defaults to 3.