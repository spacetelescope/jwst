Step Arguments
==============
The ``pathloss`` step has the following optional arguments to control
the behavior of the processing.

``--inverse`` (boolean, default=False)
  A flag to indicate whether the math operations used to apply the
  flat-field should be inverted (i.e. multiply the pathloss into
  the science data, instead of the usual division).

``--source_type`` (string, default=None)
  Force the processing to use the given source type (POINT, EXTENDED),
  instead of using the information contained in the input data. Only
  applicable to NIRSpec data.

``--user_slit_loc`` (float, default=None)
  Only applies to MIRI LRS fixed-slit exposures. Offset the target
  location along the dispersion direction of the slit by this amount,
  in units of arcsec. By definition, the center of the slit is at 0,
  and the edges in the dispersion direction are about +/-0.255 arcsec.
