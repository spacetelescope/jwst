Step Arguments
==============

The ``assign_wcs`` step has the following optional arguments to control
the behavior of the processing.

``--sip_approx`` (boolean, default=True)
  A flag to enable the computation of a SIP approximation for
  imaging modes.

``--sip_degree`` (integer, max=6, default=None)
  Polynomial degree for the forward SIP fit. "None" uses the best fit.

``--sip_max_pix_error`` (float, default=0.01)
  Maximum error for the SIP forward fit, in units of pixels. Ignored if
  ``sip_degree`` is set to an explicit value.

``--sip_inv_degree`` (integer, max=6, default=None)
  Polynomial degree for the inverse SIP fit. "None" uses the best fit.

``--sip_max_inv_pix_error`` (float, default=0.01)
  Maximum error for the SIP inverse fit, in units of pixels. Ignored if
  ``sip_inv_degree`` is set to an explicit value.

``--sip_npoints`` (integer, default=12)
  Number of points for the SIP fit.

``--slit_y_low`` (float, default=-0.55)
  Lower edge of a NIRSpec slit.

``--slit_y_high`` (float, default=0.55)
  Upper edge of a NIRSpec slit.
