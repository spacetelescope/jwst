Description
===========

Assumptions
-----------
We assume that the ``saturation`` step has already been applied to the input
science exposure, so that saturated values are appropriately flagged in the
input GROUPDQ array. We also assume that steps such as the reference pixel
correction (``refpix``) and non-linearity correction (``linearity``) have been applied, so
that the input data ramps do not have any non-linearities or noise above the modeled Poission
and read noise due to instrumental effects. The absence of any of these preceding corrections
or residual non-linearities or noise can lead to the false detection of jumps in the ramps,
due to departure from linearity.

The ``jump`` step will automatically skip execution if the input data contain fewer
than 5 groups per integration, because the baseline algorthim requires four first
differences to work.

Algorithm
---------
This routine detects jumps in an exposure by looking for outliers
in the up-the-ramp signal for each pixel in each integration within
an input exposure. On output, the GROUPDQ array is updated with the DQ flag
"JUMP_DET" to indicate the location of each jump that was found.
In addition, any pixels that have non-positive or NaN values in the gain
reference file will have DQ flags "NO_GAIN_VALUE" and "DO_NOT_USE" set in the
output PIXELDQ array.
The SCI and ERR arrays of the input data are not modified.

The current implementation uses the two-point difference method described
in Anderson&Gordon2011_.

Two-Point Difference Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The two-point difference method is applied to each integration as follows:

* Compute the first differences for each pixel (the difference between
  adjacent groups)
* Compute the clipped (dropping the largest difference) median of the first differences for each pixel.
* Use the median to estimate the Poisson noise for each group and combine it
  with the read noise to arrive at an estimate of the total expected noise for
  each difference.
* Compute the "difference ratio" as the difference between the first differences
  of each group and the median, divided by the expected noise.
* If the largest "difference ratio" is greater than the rejection threshold,
  flag the group corresponding to that ratio as having a jump.
* If a jump is found in a given pixel, iterate the above steps with the
  jump-impacted group excluded, looking for additional lower-level jumps
  that still exceed the rejection threshold.
* Stop iterating on a given pixel when no new jumps are found

Note that any ramp values flagged as SATURATED in the input GROUPDQ array
are not used in any of the above calculations and hence will never be
marked as containing a jump.

Multiprocessing
===============
This step has the option of running in multiprocessing mode. In that mode it will
split the input data cube into a number of row slices based on the number of available
cores on the host computer and the value of the max_cores input parameter. By
default the step runs on a single processor. At the other extreme if max_cores is
set to 'all', it will use all available cores (real and virtual). Testing has shown
a reduction in the elapsed time for the step proportional to the number of real
cores used. Using the virtual cores also reduces the elasped time but at a slightly
lower rate than the real cores.

If multiprocessing is requested the input cube will be divided into a number of
slices in the row dimension (with the last slice being slightly larger, if needed).
The slices are then sent to twopoint_difference.py by detect_jumps. After all the
slices have finished processing, detect_jumps assembles the output group_dq cube
from the slices.

Subarrays
=========
The use of the reference files is flexible. Full-frame reference
files can be used for all science exposures, in which case subarrays will be
extracted from the reference file data to match the science exposure, or
subarray-specific reference files may be used.

.. _Anderson&Gordon2011: https://ui.adsabs.harvard.edu/abs/2011PASP..123.1237A
