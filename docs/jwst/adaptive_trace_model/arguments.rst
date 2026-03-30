Step Arguments
==============

The ``adaptive_trace_model`` step has the following step-specific arguments:

``--fit_threshold`` (float, default=10.0)
  Signal threshold value for attempting a spline fit, specified as the number of
  standard deviations from the overall mean for the image. This parameter controls
  whether a spline fit is attempted for a given spectral region (slit or slice).
  Higher values will create spline models for fewer slices; lower values will attempt
  to fit more slices. If set to 0, all slices will be fit.

``--oversample`` (float, default=1.0)
  Use the trace model to oversample the data by this factor.  If the value is 1.0,
  the only change to the datamodel is to attach a spectral trace image, in the
  ``trace_model`` attribute.  If any other value is specified, the input flux, error,
  and variance images are replaced with interpolated data oversampled onto a new pixel
  grid.

``--slope_limit`` (float, default=0.1)
  Slope limit for using the spline model to compute the trace.  This parameter is used to
  distinguish between bright, compact sources for which the spline model is appropriate and
  faint, diffuse sources for which it is not.
  For compact sources (high slope), the spline models are saved in the output trace model
  and are used to compute the oversampled flux if oversampling is performed.
  For diffuse sources (low slope), a linear interpolation is used in oversampling.
  Set the slope limit to lower values to use the spline model for fainter sources.
  If set to zero, the spline model will always be used.

``--psf_optimal`` boolean(default=False)
  If set to True, the values for ``fit_threshold`` and ``slope_limit`` are ignored and
  the spline models are fit and used for all data.  Also, residual differences from
  the spline model are not interpolated and added to the spline fits in oversampling.
  This option is generally only appropriate for simple, isolated point sources.

``--save_intermediate_results`` boolean(default=False)
  If set to True, additional image models are saved to disk for inspection, containing the full
  spline model, the spline model as used for compact sources, the linearly interpolated
  data, and the residual flux model.  The saved files will have suffix
  ``spline_full``, ``spline_used``, ``linear_interp``, and ``spline_residual``,
  respectively.  The linear interpolation and residual flux files are saved only
  if oversampling was performed.

``--maximum_cores``: The number of available cores that will be
  used for multi-processing in this step. The default value is '1', which does not use
  multi-processing. The other options are either an integer, 'quarter', 'half', or 'all'.
  Note that these fractions refer to the total available cores and on most CPUs these include
  physical and virtual cores.  Note that only the spline fitting portion of the code will
  use multiprocessing. The oversampling portion is not affected by this flag. Thus the
  speed up will be less than linear with the number of cores, since only a portion of
  the code is affected. The relative fractional speedup will be largest for cases in
  which the fit_threshold and slope_limit parameter values, and/or the psf_optimal flag,
  result in a larger fraction of the image having the spline fit performed.
