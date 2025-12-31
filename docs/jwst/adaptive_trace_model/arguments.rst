Step Arguments
==============

The ``adaptive_trace_model`` step has the following step-specific arguments:

``--fit_threshold`` (float, default=10.0)
  Limiting sigma for fitting splines. This parameter controls whether a spline fit is
  attempted at all for a spectral region (slit or slice).  Higher values will create
  spline models for fewer slices; lower values will attempt to fit more slices.

``--oversample`` (float, default=1.0)
  Use the trace model to oversample the data by this factor.  If the value is 1.0,
  the only change to the datamodel is to attach a spectral trace image, in the
  ``trace_model`` attribute.  If any other value is specified, the input flux, error,
  and variance images are replaced with interpolated data oversampled onto a new pixel
  grid.

``--slope_limit`` (float, default=0.1)
  Slope limit for using splines in oversample.  This parameter is used to distinguish
  between bright, compact sources and faint diffuse sources for oversampling purposes.
  For compact sources (high slope), the spline models are used in the interpolation.
  For diffuse sources (low slope), a linear interpolation is used.  Set the slope
  limit to lower values to use the spline model for fainter sources.
