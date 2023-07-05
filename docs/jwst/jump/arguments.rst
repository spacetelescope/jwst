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

  After a jump of at least 'after_jump_flag_dn1' DN, groups up to 'after_jump_flag_time1'
  seconds will be also flagged as jumps. That pair of arguments is defined as:
* ``--after_jump_flag_dn1``: A floating point value in units of DN
* ``--after_jump_flag_time1``: A floating point value in units of seconds

  A second threshold and time can also be set: after a jump of at least 'after_jump_flag_dn2' DN,
  groups up to 'after_jump_flag_time2' seconds will be also flagged as jumps. That pair of arguments
  is defined as:
* ``--after_jump_flag_dn2``: A floating point value in units of DN
* ``--after_jump_flag_time2``: A floating point value in units of seconds

* ``--expand_large_events``:  A boolean parameter that controls whether the jump step will expand the number of pixels that are flagged around large cosmic ray events. These are know as "snowballs" in the near-infrared detectors and "showers" for the MIRI detectors. In general, this should be set to True.

* ``--min_jump_area``: The minimum number of contiguous pixels needed to trigger the expanded flagging of large cosmic rays events.

* ``--min_sat_area``:  The minimum number of saturated pixels required to meet "sat_required_snowball".

* ``--expand_factor``: A multiplicative factor applied to the enclosing ellipse for snowballs. This larger area will have all pixels flagged as having a jump.

* ``--use_ellipses``:  deprecated

* ``--sat_required_snowball``: A boolean value that if True requires that there are saturated pixels within the enclosed jump circle.

* ``--min_sat_radius_extend``: The minimum radius of the saturated core of a snowball required to for the radius of the saturated core to be extended.

* ``--sat_expand``: Number of pixels to add to the radius of the saturated core of snowballs

* ``--edge_size``: The distance from the edge of the detector where saturated cores are not required for snowball detection

* ``--find_showers``: Turn on the detection of showers for the MIRI detectors

* ``--extend_snr_threshold``: The SNR minimum for the detection of faint extended showers in MIRI

* ``--extend_min_area``: The required minimum area of extended emission after convolution for the detection of showers in MIRI

* ``--extend_inner_radius``: The inner radius of the ring_2D_kernel that is used for the detection of extended emission in showers

* ``--extend_outer_radius``: The outer radius of the Ring2DKernal that is used for the detection of extended emission in showers

* ``--extend_ellipse_expand_ratio``: Multiplicative factor to expand the radius of the ellipse fit to the detected extended emission in MIRI showers

* ``--time_masked_after_showers``: Number of seconds to flag groups as jump after a detected extended emission in MIRI showers

* ``--max_extended_radius``: The maxiumum extension of the jump and saturation that will be flagged for showers or snowballs

* ``--minimum_groups``: The minimum number of groups to run the jump step

* ``--minimum_sigclip_groups``: The minimum number of groups to switch the jump detection to use sigma clipping

* ``--only_use_ints``: If true the sigma clipping is applied only for a given group across all ints. If not, all groups from all ints are used for the sigma clipping.


