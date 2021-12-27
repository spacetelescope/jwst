Step Arguments
==============

The ``flat_field`` step has the following optional arguments to control
the behavior of the processing.

``--save_interpolated_flat`` (boolean, default=False)
  A flag to indicate whether to save to a file the NIRSpec
  flat field that was constructed on-the-fly by the step.
  Only relevant for NIRSpec data.

``--user_supplied_flat`` (string, default=None)
  The name of a user-supplied flat-field reference file.

``--inverse`` (boolean, default=False)
  A flag to indicate whether the math operations used to apply the
  flat-field should be inverted (i.e. multiply the flat-field into
  the science data, instead of the usual division).
