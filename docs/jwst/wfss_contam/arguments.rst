.. _wfss_contam_step_args:

Step Arguments
==============
The ``wfss_contam`` step uses the following optional arguments.

``--save_simulated_image`` (bool, default=False)
  Whether the full-frame simulated grism image containing all
  simulated spectra within the field-of-view should be saved to a file. The file
  name uses a product type suffix of "simul".

``--save_contam_images`` (bool, default=False)
  This parameter is deprecated and has no effect.

``--maximum_cores`` (str, default='1')
  Number of available cores that will be
  used for multi-processing in this step. The default value is '1', which does not use
  multi-processing. The other options are either an integer, 'quarter', 'half', or 'all'.
  Note that these fractions refer to the total available cores and on most CPUs these include
  physical and virtual cores.

``--orders`` (list, default=None)
  Which grism orders to simulate. The default value is None, which
  means all orders that are defined in the wavelength and specwcs reference files
  for that instrument will be simulated.
  To specify a single order from the command line, use e.g. "0," or "1,"
  (the comma allows the code to identify this as a list).

``--magnitude_limit`` (float, default=None)
  Magnitude limit for sources to be included in the simulation, taken from the input source
  catalog's isophotal AB magnitude column.
  The limit is scaled according to the relative sensitivity of each spectral order, such that
  fewer sources are included in orders with lower sensitivity.
  The default value of None means no magnitude limit is applied to any of the orders
  and all sources are included.

``--wl_oversample`` (int, default=2)
  Indicates the oversampling factor for the wavelength grid used in the
  simulation of the dispersed spectra.

``--max_pixels_per_chunk`` (int, default=5000)
  Sets the maximum number of direct image pixels to run through the grism transforms at once.
  Decreasing this value will typically reduce the memory usage of the step. The effect on runtime
  depends on the machine hardware and whether multi-processing is enabled.

Polynomial fitting parameters
-----------------------------

``--polyfit_degree`` (int, default=None)
  Maximum degree of the polynomial used to fit the spectral shape of each
  source when computing the contamination estimate. If ``None``, no polynomial
  fitting is applied and the direct image flux values are used as-is (i.e., a "flat" spectrum).

``--n_iterations`` (int, default=1)
  Number of polynomial fitting iterations to perform. On each
  iteration the contamination estimate is recomputed using the flux-modeled spectra from the
  previous iteration. The default of 1 means a single, non-iterative correction. Has no effect if
  ``polyfit_degree`` is ``None``.

``--l2_alpha`` (float, default=0.1)
  Alpha parameter for L2 regularization in the polynomial fitting.
  This is used to prevent the fitted coefficients from blowing up in cases of severe contamination.

``--rejection_threshold`` (float, default=0.1)
  Threshold for rejecting polynomial fits based on the fitted constant term coefficient.
  If the absolute value of that coefficient deviates from unity by more than this threshold, the fit is
  rejected and the contamination estimate for that source is not updated on that iteration.
  This is used to avoid fits "blowing up" in cases where the polynomial fit has returned
  an unphysical total flux level, which typically occurs if background subtraction was imperfect 
  or if the source sits in a highly contaminated region. If set to None, no rejection will be performed.
