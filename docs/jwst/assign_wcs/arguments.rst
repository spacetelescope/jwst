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
  Lower edge of a NIRSpec slit in slit units (see below).

``--slit_y_high`` (float, default=0.55)
  Upper edge of a NIRSpec slit in slit units.

  Slit units are specified as a fraction of the nominal slit or shutter height.
  In these relative units, -0.5 is the nominal *top* edge of the slit open area in
  the image frame, 0.0 is the center of the slit, and 0.5 is the nominal *bottom* edge
  of the slit open area in the image frame.

  Set the `slit_y_low` value to a larger negative value to include more pixels
  at the *top* of a slit bounding box; a smaller negative value to to include fewer
  pixels.

  Set the `slit_y_high` value to a larger positive value to include more pixels
  at the *bottom* of a slit bounding box; a smaller positive value to to include fewer
  pixels.

  For MSA slits in NIRSpec MOS mode, the slit units are relative to a single open shutter
  height, not a full "slitlet" composed of multiple shutters.  For this mode,
  the `slit_y_low` and `slit_y_high` values determine the padding of the slitlet edges relative
  to the centers of the edge shutters. An additional margin of half a shutter is also
  automatically added to the slit boundaries specified by these parameters, corresponding
  to a couple extra detector pixels at the top and bottom of each slitlet.

  Please note that significantly expanding the slit limit values may introduce
  contamination from adjacent open slits when the slit images are extracted
  in the :ref:`extract_2d <extract_2d_step>` step, and/or may expand the slits
  into regions that cannot be flat fielded in the :ref:`flat_field <flatfield_step>`
  step.  These parameters should be used with caution.

``--nrs_ifu_slice_wcs`` (boolean, default=False)

  If True and the exposure type is NIRSpec IFU, then a full slice-based
  WCS that propagates slice IDs is produced.  This is intended primarily for
  diagnostic purposes.  If False and the exposure type is NIRSpec IFU,
  a slice map is internally applied to produce a fully coordinate-based
  WCS pipeline that does not require slice IDs on input.
