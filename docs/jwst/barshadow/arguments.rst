Step Arguments
==============
The ``barshadow`` step has the following optional arguments.

``--inverse`` (boolean, default=False)
  A flag to indicate whether the math operations used to apply the
  correction should be inverted (i.e. multiply the correction into
  the science data, instead of the usual division).

``--source_type`` (string, default=None)
  Force the processing to use the given source type (POINT, EXTENDED),
  instead of using the information contained in the input data.
