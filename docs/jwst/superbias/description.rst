Description
===========

:Class: `jwst.superbias.SuperBiasStep`
:Alias: superbias

The superbias subtraction step removes the fixed detector bias from a
science data set by subtracting a superbias reference image.

Algorithm
---------
The 2-D superbias reference image is subtracted from every group in every
integration of the input science ramp data. Any NaN's that might be present
in the superbias image are set to a value of zero before being subtracted
from the science data, such that those pixels effectively receive no correction.

The DQ array from the superbias reference file is combined with the science
exposure "PIXELDQ" array using a bitwise OR operation.

The ERR and GROUPDQ arrays in the science ramp data are unchanged.

NIRCam Frame 0
--------------
If the NIRCam frame zero data cube is present in the input data, the
image for each integration has the superbias reference image subtracted
from it, in the same way as the regular science data.

Subarrays
---------
If the subarray mode of the superbias reference file matches that of the
science exposure, the reference data are directly subtracted. If the superbias
reference file contains full-frame data, while the science exposure is a
subarray mode, the corresponding subarray is extracted from the superbias
reference data before being subtracted.

NIRSpec IRS2
------------
No special handling is necessary for science data taken in the IRS2 readout
mode, because matching IRS2 superbias reference files are supplied in CRDS.
