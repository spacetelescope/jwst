Description
============

:Class: `jwst.linearity.LinearityStep`
:Alias: linearity

Assumptions
-----------
It is assumed that the input science exposure data for near-IR instruments
have had the :ref:`superbias <superbias_step>` subtraction step applied,
therefore the correction coefficients stored in the linearity reference files
for those instruments must have been derived from data that has also been
bias subtracted.
MIRI data, on the other hand, do not receive bias subtraction
(see :ref:`calwebb_detector1 <calwebb_detector1>`) and hence the linearity
correction is derived from data that has not been bias subtracted.

It is also assumed that the saturation step has already been applied to
the input data, so that saturation flags are set in the GROUPDQ array of
the input science data.

Algorithm
---------
The algorithm for this step is called from the external package ``stcal``, an STScI
effort to unify common calibration processing algorithms for use by multiple observatories.

The linearity step applies the "classic" linearity correction adapted from
the HST WFC3/IR linearity correction routine, correcting science data values
for detector non-linearity. The correction is applied pixel-by-pixel,
group-by-group, integration-by-integration within a science exposure.

The correction is represented by an nth-order polynomial for
each pixel in the detector, with n+1 arrays of coefficients read from the
linearity reference file.

The algorithm for correcting the observed pixel value in each group of an
integration is currently of the form:

.. math::
   F_\text{c} = c_{0} + c_{1}F + c_{2}F^2 + c_{3}F^3 + ... + c_{n}F^n

where :math:`F` is the observed counts (in DN), :math:`c_n` are the polynomial
coefficients, and :math:`F_\text{c}` is the corrected counts. There is no
limit to the order of the polynomial correction; all coefficients contained in
the reference file will be applied.

Upon successful completion of the linearity correction the S_LINEAR keyword is
set to "COMPLETE".

Special Handling
++++++++++++++++

- Pixels having at least one correction coefficient equal to NaN will not have
  the linearity correction applied and the DQ flag "NO_LIN_CORR" is added to
  the science exposure PIXELDQ array.

- Pixels that have the "NO_LIN_CORR" flag set in the DQ array of the linearity
  reference file will not have the correction applied and the "NO_LIN_CORR" flag
  is added to the science exposure PIXELDQ array.

- Pixel values that have the "SATURATED" flag set in a particular group of the
  science exposure GROUPDQ array will not have the linearity correction
  applied to that group. Any groups for that pixel that are not flagged as
  saturated will be corrected.

The ERR array of the input science exposure is not modified.

The flags from the linearity reference file DQ array are propagated into the
PIXELDQ array of the science exposure using a bitwise OR operation.

NIRCam Frame 0
--------------

If the NIRCam "Frame 0" data are included in the input, the linearity
correction is applied to each integration's frame zero image in the same
way as it's applied to the normal science data cube. The corrected
frame zero data are returned as part of the overall datamodel being
processed.

Subarrays
---------

This step handles input science exposures that were taken in subarray modes
in a flexible way. If the reference data arrays are the same size as the
science data, they will be applied directly. If there is a mismatch, the
routine will extract a matching subarray from the reference file data arrays
and apply them to the science data. Hence full-frame reference files can be
used for both full-frame and subarray science exposures, or
subarray-dependent reference files can be provided if desired.
