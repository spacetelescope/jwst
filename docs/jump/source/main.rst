Description
===========

Assumptions
-----------
We assume that the saturation step has already been applied to the input
science exposure, so that saturated values are appropriately flagged in the
input GROUPDQ array. We also assume that steps such as the reference pixel
correction (bias drift) and non-linearity correction have been applied, so
that the input data ramps do not have any non-linearities due to instrumental
effects. The absence of any of these preceding corrections can lead to the
false detection of jumps in the ramps, due to departure from linearity.

The step will automatically skip execution if the input data contain fewer
than 3 groups per integration, because it's impossible to detect jumps with
only 1 or 2 groups.

Algorithms
----------
This routine detects jumps in an exposure by looking for outliers
in the up-the-ramp signal for each pixel in each integration within
an input exposure. On output, the GROUPDQ array of the data is updated to
reflect the location of each jump that was found. The SCI, ERR, and PIXELDQ
arrays of the input data are not modified.

The current implementation uses a combination of the two-point difference
and y-intercept methods described in Anderson and Gordon, PASP 132, 1237
(2011). These methods are applied to the data in two passes. The
two-point difference method, which is computationally cheap (fast), is
applied in a first pass to all pixels in each integration. The y-intercept
method, which is computationally expensive (really slow), is applied in a
second pass to only those pixels that have signal levels in the read noise
regime, where it provides somewhat better results than the two-point
difference method. The y-intercept method is only applied if the jump step
parameter 'do_yintercept' is set to True (default is False).

Two-Point Difference Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The two-point difference method is applied to each integration as follows:

* Compute the first differences for each pixel (the difference between
  adjacent groups)
* Compute the median of the first differences for each pixel
* Use the median to estimate the Poisson noise for each group and combine it
  with the read noise to arrive at an estimate of the total expected noise for
  each group
* Take the ratio of the first differences and the total noise for each group
* If the largest ratio is above the rejection threshold, flag the group
  corresponding to that ratio as having a jump
* If a jump is found, iterate the above steps with the jump-impacted group
  excluded, looking for additional jumps
* Stop iterating on a given pixel when no new jumps are found

Y-Intercept Method
^^^^^^^^^^^^^^^^^^
After the two-point difference method has been applied, the slope of each
pixel is calculated, taking into account any jumps that were flagged. The
y-intercept method is then applied, if requested, to all pixels with 
signals (e/sec) that fall below the value of the ``yint_threshold`` parameter
value.

The y-intercept method is applied as follows:

* For each group within a pixel ramp compute the slope and y-intercept 
  for two portions of the ramp: all groups preceding the current group and
  all groups following the current group. The slope and y-intercept are
  computed using a full covariance matrix scheme that takes into account the
  correlated Possion noise and uncorrelated read noise for each group.
* Compute the difference between the y-intercepts for each of the
  above pairs.
* Compute the expected uncertainty in the y-intercepts based on the Poisson
  noise and read noise in each group.
* Take the ratio of the y-intercept differences and the expected uncertainties.
* If the largest of these ratios for a given pixel is above the rejection
  threshold, flag the group corresponding to that ratio as having a jump
* If a jump is found, split the ramp into two segments (one preceding and one
  following the jump) and iterate the above steps on each resulting
  semi-ramp, looking for additional jumps
* Stop iterating on a given pixel when no new jumps are found

Step Arguments
==============
The Jump step has three optional arguments that can be set by the user:

* ``--rejection_threshold``: A floating-point value that sets the sigma
  threshold for jump detection.
* ``--do_yintercept``: A True/False value that indicates whether to apply
  the y-intercept method.
* ``--yint_threshold``: A floating-point value that sets the signal
  threshold for applying the y-intercept method to individual pixels.

Subarrays
---------
The use of the reference files is flexible. Full-frame reference
files can be used for all science exposures, in which case subarrays will be
extracted from the reference file data to match the science exposure, or
subarray-specific reference files may be used.
