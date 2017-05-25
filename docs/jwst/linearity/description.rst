
Description
============

Assumptions
-----------

Beginning with the Build 5 pipeline, it is assumed that the input science
exposure data from near-IR instruments have had the superbias subtraction
applied, therefore the correction coefficients stored in the linearity
reference files for those instruments must also have been
derived from data that had the zero group subtracted.

It is also assumed that the saturation step has already been applied to
the science data, so that saturation flags are set in the GROUPDQ array of
the input science data.

Algorithm
---------

The linearity step applies the "classic" linearity correction adapted from
the HST WFC3/IR linearity correction routine, correcting science data values
for detector non-linearity. The correction is applied pixel-by-pixel,
group-by-group, integration-by-integration within a science exposure.  Pixels
having at least one correction coefficient equal to NaN (not a number), or are
flagged with "Linearity Correction not determined for pixel" (NO_LIN_CORR) in
the PIXELDQ
array will not have the linearity correction applied. Pixel values flagged as
saturated in the GROUPDQ array for a given group will not have the linearity
correction applied. All non-saturated groups for such a pixel will have the
correction applied.

The correction is represented by an nth-order polynomial for
each pixel in the detector, with n+1 arrays of coefficients read from the
linearity reference file.

The algorithm for correcting the observed pixel value in each group of an
integration is currently of the form
:math:`F_\text{c} = c_{0} + c_{1}F + c_{2}F^2 + c_{3}F^3 ...`

where :math:`F` is the observed counts (in DN), :math:`c_n` are the polynomial
coefficients, and :math:`F_\text{c}` is the corrected counts. There is no
limit to the order of the polynomial correction; all coefficients contained in
the reference file will be applied.

The ERR array of the input science exposure is not modified.

The values from the linearity reference file DQ array are propagated into the
PIXELDQ array of the input science exposure using a bitwise OR operation.

Subarrays
---------

This step handles input science exposures that were taken in subarray modes
in a flexible way. If the reference data arrays are the same size as the
science data, they will be applied directly. If there is a mismatch, the
routine will extract a matching subarray from the reference file data arrays
and apply them to the science data. Hence full-frame reference files can be
used for both full-frame and subarray science exposures, or
subarray-dependent reference files can be provided if necessary.
