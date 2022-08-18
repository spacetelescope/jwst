Arguments
=========

The ``jump`` step has five optional arguments that can be set by the user:

* ``--rejection_threshold``: A floating-point value that sets the sigma
  threshold for jump detection. In the code sigma is determined using the read noise from the
  read noise reference file and the Poisson noise (based on the median difference between
  samples, and the gain reference file). Note that any noise source beyond these two that
  may be present in the data will lead to an increase in the false positive rate and thus
  may require an increase in the value of this parameter. The default value of 4.0 for the
  rejection threshold will yield 6200 false positives for every million pixels, if the noise
  model is correct.

* ``--maximum_cores``: The fraction of available cores that will be
  used for multi-processing in this step. The default value is 'none' which does not use
  multi-processing. The other options are 'quarter', 'half', and 'all'. Note that these
  fractions refer to the total available cores and on most CPUs these include physical
  and virtual cores. The clock time for the step is reduced
  almost linearly by the number of physical cores used on all machines. For example, on an Intel CPU with
  six real cores and 6 virtual cores setting maximum_cores to 'half' results in a
  decrease of a factor of six in the clock time for the step to run. Depending on the system
  the clock time can also decrease even more with maximum_cores is set to 'all'.

* ``--flag_4_neighbors``: If set to True (default is True) it will cause the four perpendicular
  neighbors of all detected jumps to be flagged as a jump. This is needed because of
  the inter-pixel capacitance (IPC) causing a small jump in the neighbors. The small jump
  might be below the rejection threshold but will affect the slope determination of
  the pixel. The step will take about 40% longer to run when this is set to True.

* ``--max_jump_to_flag_neighbors``: A floating point value in units of sigma that limits
  the flagging of neighbors. Any jump above this cutoff will not have its neighbors flagged.
  The concept is that the jumps in neighbors will be above the rejection-threshold and thus
  be flagged as primary jumps. The default value is 200.

* ``--min_jump_to_flag_neighbors``: A floating point value in units of sigma that limits
  the flagging of neighbors of marginal detections. Any primary jump below this value will
  not have its neighbors flagged. The goal is to prevent flagging jumps that would be too
  small to significantly affect the slope determination.  The default value is 10.

*  ``--after_jump_flag_dn1``: A floating point value in units of DN that gives the
   minimum ramp jump for after jump flagging for the 1st after ramp flagging.

* ``--after_jump_flag_time1``: A floating point value in units of seconds that gives the
  time since the ramp jump to flag groups as ramp jumps for the 1st after jump flagging.

*  ``--after_jump_flag_dn1``: A floating point value in units of DN that gives the
   minimum ramp jump for after jump flagging for the 2nd after ramp flagging.

* ``--after_jump_flag_time1``: A floating point value in units of seconds that gives the
  time since the ramp jump to flag groups as ramp jumps for the 2nd after jump flagging.
