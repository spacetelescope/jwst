Description
===========

Assumptions
-----------
We assume that the ``saturation`` step has already been applied to the input
science exposure, so that saturated values are appropriately flagged in the
input GROUPDQ array. We also assume that steps such as the reference pixel
correction (``refpix``) and non-linearity correction (``linearity``) have been applied, so
that the input data ramps do not have any non-linearities due to instrumental
effects. The absence of any of these preceding corrections can lead to the
false detection of jumps in the ramps, due to departure from linearity.

The ``jump`` step will automatically skip execution if the input data contain fewer
than 3 groups per integration, because it's impossible to detect jumps with
only 1 or 2 groups.

Algorithm
---------
This routine detects jumps in an exposure by looking for outliers
in the up-the-ramp signal for each pixel in each integration within
an input exposure. On output, the GROUPDQ array is updated with the DQ flag
JUMP_DET to indicate the location of each jump that was found.
In addition, any pixels that have non-positive or NaN values in the gain
reference file will have DQ flags NO_GAIN_VALUE and DO_NOT_USE set in the
output PIXELDQ array.
The SCI and ERR arrays of the input data are not modified.

The current implementation uses the two-point difference method described
in Anderson and Gordon, PASP 132, 1237 (2011).


Two-Point Difference Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The two-point difference method is applied to each integration as follows:

* Compute the first differences for each pixel (the difference between
  adjacent groups)
* Compute the clipped median of the first differences for each pixel
* Use the median to estimate the Poisson noise for each group and combine it
  with the read noise to arrive at an estimate of the total expected noise for
  each group
* Compute the "difference ratio" as the difference between the first differences
  of each group and the median, divided by the expected noise
* If the largest "difference ratio" is greater than the rejection threshold,
  flag the group corresponding to that ratio as having a jump
* If a jump is found in a given pixel, iterate the above steps with the
  jump-impacted group excluded, looking for additional lower-level jumps
  that still exceed the rejection threshold
* Stop iterating on a given pixel when no new jumps are found

Note that any ramp values flagged as SATURATED in the input GROUPDQ array
are not used in any of the above calculations and hence will never be
marked as containing a jump.

Step Arguments
==============
The ``jump`` step has one optional argument that can be set by the user:

* ``--rejection_threshold``: A floating-point value that sets the sigma
  threshold for jump detection.


Subarrays
=========
The use of the reference files is flexible. Full-frame reference
files can be used for all science exposures, in which case subarrays will be
extracted from the reference file data to match the science exposure, or
subarray-specific reference files may be used.
