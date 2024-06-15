Arguments
=========
The ramp fitting step has the following optional arguments that can be set by the user:

* ``--algorithm``: A string to select the desired algorithm.  The available
  values are "OLS" to select the python implementation of the Ordinary
  Least Squares algorithm and  "OLS_C" to select the C extension
  implementation of OLS.  The algorithm defaults to "OLS".

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

* ``--maximum_cores``: The number of available cores that will be
  used for multi-processing in this step. The default value is '1', which results in no
  multi-processing. Other options are either an integer, 'quarter', 'half', and 'all'.
  Note that these fractions refer to the total available cores and on most CPUs these include
  physical and virtual cores. The clock time for the step is reduced almost linearly by the
  number of physical cores used on all machines. For example, on an Intel CPU with
  six real cores and six virtual cores, setting maximum_cores to 'half' results in a
  decrease of a factor of six in the clock time for the step to run. Depending on the system,
  the clock time can also decrease even more with maximum_cores set to 'all'.
  Setting the number of cores to an integer can be useful when running on machines with a
  large number of cores where the user is limited in how many cores they can use.
