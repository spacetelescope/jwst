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

Algorithm
---------
This routine detects jumps in an exposure by looking for outliers
in the up-the-ramp signal for each pixel in each integration within
an input exposure. On output, the GROUPDQ array of the data is updated to
reflect the location of each jump that was found. The SCI, ERR, and PIXELDQ
arrays of the input data are not modified.

The current implementation uses the two-point difference method described 
in Anderson and Gordon, PASP 132, 1237 (2011). 


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


Step Arguments
==============
The Jump step has one optional argument that can be set by the user:

* ``--rejection_threshold``: A floating-point value that sets the sigma
  threshold for jump detection.


Subarrays
---------
The use of the reference files is flexible. Full-frame reference
files can be used for all science exposures, in which case subarrays will be
extracted from the reference file data to match the science exposure, or
subarray-specific reference files may be used.
