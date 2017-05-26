
Description
===========

The superbias subtraction step removes the fixed detector bias from a
science data set by subtracting a superbias reference image.

Algorithm
---------

The 2-D superbias reference image is subtracted from every group in every
integration of the input science ramp data. Any NaN's that might be present
in the superbias image are set to a value of zero before being subtracted
from the science data, such that those pixels effectively receive no correction.

The DQ array from the superbias reference file is combined with the science
exposure PIXELDQ array using a bit-wise OR operation.

The ERR arrays in the science ramp data are unchanged.

Subarrays
---------

If the subarray mode of the superbias reference file matches that of the
science exposure, the reference data are applied as-is. If the superbias
reference file contains full-frame data, while the science exposure is a
subarray mode, the corresponding subarray will be extracted from the superbias
reference data before being applied.
