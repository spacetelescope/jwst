Description
============

The ``saturation`` step flags pixels at or below the A/D floor or above the
saturation threshold.  Pixels values are flagged as saturated if the pixel value
is larger than the defined saturation threshold.  Pixel values are flagged as
below the A/D floor if they have a value of zero DN.

This step loops over all integrations within an exposure, examining each one
group-by-group, comparing the pixel values in the SCI array with defined
saturation thresholds for each pixel. When it finds a pixel value in a given
group that is above the saturation threshold (high saturation), it sets the
"SATURATED" flag in the corresponding location of the "GROUPDQ" array in the
science exposure.  When it finds a pixel in a given group that has a zero or
negative value (below  the A/D floor), it sets the "AD_FLOOR" and "DO_NOT_USE"
flags in the corresponding location of the "GROUPDQ" array in the science
exposure  For the saturation case, it also flags all subsequent groups for that
pixel as saturated. For example, if there are 10 groups in an integration and
group 7 is the first one to cross the saturation threshold for a given pixel,
then groups 7 through 10 will all be flagged for that pixel.

Pixels with thresholds set to NaN or flagged as "NO_SAT_CHECK" in the saturation
reference file have their thresholds set to the 16-bit A-to-D converter limit
of 65535 and hence will only be flagged as saturated if the pixel reaches that
hard limit in at least one group. The "NO_SAT_CHECK" flag is propagated to the
PIXELDQ array in the output science data to indicate which pixels fall into
this category.

NIRSpec data acquired using the "IRS2" readout pattern require special
handling in this step, due to the extra reference pixel values that are interleaved
within the science data. The saturation reference file data does not contain
extra entries for these pixels. The step-by-step process is as follows:

- Retrieve and load data from the appropriate "SATURATION" reference file from CRDS

- If the input science exposure used the NIRSpec IRS2 readout pattern:

 * Create a temporary saturation array that is the same size as the IRS2 readout

 * Copy the saturation threshold values from the original reference data into
   the larger saturation array, skipping over the interleaved reference pixel
   locations within the array

- If the input science exposure used a subarray readout, extract the matching
  subarray from the full-frame saturation reference file data

- For pixels that contain NaN in the reference file saturation threshold array
  or are flagged in the reference file with "NO_SAT_CHECK" (no saturation check
  available), propagate the "NO_SAT_CHECK" flag to the science data PIXELDQ array

- For each group in the input science data, set the "SATURATION" flag in the
  "GROUPDQ" array if the pixel value is greater than or equal to the saturation
  threshold from the reference file

Subarrays
=========
The ``saturation`` step will accept either full-frame or subarray saturation reference files.
If only a full-frame reference file is available, the step will extract a
subarray to match that of the science exposure. Otherwise, subarray-specific
saturation reference files will be used if they are available.
