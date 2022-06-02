Arguments
=========
The ramp fitting step has three optional arguments that can be set by the user:

* ``--save_opt``: A True/False value that specifies whether to write
  the optional output product. Default is False.

* ``--opt_name``: A string that can be used to override the default name
  for the optional output product.

* ``--int_name``: A string that can be used to override the default name
  for the per-integration product.

* ``--suppress_one_group``: A boolean to suppress computations for saturated ramps
  with only one good (unsaturated) sample.  The default is set to True to suppress these computations,
  which will compute all values for the ramp the same as if the entire ramp were
  saturated.

* ``--maximum_cores``: The fraction of available cores that will be
  used for multi-processing in this step. The default value is 'none' which does not use
  multi-processing. The other options are 'quarter', 'half', and 'all'. Note that these
  fractions refer to the total available cores and on most CPUs these include physical
  and virtual cores. The clock time for the step is reduced
  almost linearly by the number of physical cores used on all machines. For example, on an Intel CPU with
  six real cores and 6 virtual cores setting maximum_cores to 'half' results in a
  decrease of a factor of six in the clock time for the step to run. Depending on the system
  the clock time can also decrease even more with maximum_cores is set to 'all'.
