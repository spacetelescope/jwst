Description
===========

:Class: `jwst.ramp_fitting.RampFitStep`
:Alias: ramp_fit

This step determines the mean count rate, in units of counts per second, for
each pixel by performing a linear fit to the "up the ramp" data in the input file.
All fitting is done using the "ordinary least squares" (OLS) method.
The fit is performed independently for each pixel. Bad values flagged via
DQ flags are rejected from the fits.

The input to the step is assumed to be the fully-corrected and flagged 4-D
data resulting from applying all previous steps of the
:ref:`calwebb_detector1 <calwebb_detector1>` pipeline and will nominally be
the output from the :ref:`jump detection <jump_step>` step. It is in fact
vital that all anomalies such as saturation, non-linearity, and CR jumps
be corrected or appropriately flagged in order to obtain useful results
from ramp fitting.

There are two output products created by default, with a third optional
product also available:

#. The primary output file ("rate") contains slope and variance/error
   estimates for each pixel that are the result of averaging over all
   integrations in the exposure. This is a product with 2-D data arrays.
#. The secondary product ("rateints") contains slope and variance/error
   estimates for each pixel on a per-integration basis, stored as 3-D
   data cubes.
#. The third, optional, output product contains detailed
   fit information for every ramp segment for each pixel.

The three types of output products are described in more detail below.

The count rate for each pixel is determined by a linear fit to the
cosmic-ray-free and saturation-free ramp intervals for each pixel, with any
intervening groups flagged as "DO_NOT_USE" excluded from the fits. Hereafter
such intervals will be referred to as a ramp "segment." The fitting algorithm uses an 
'optimal' weighting scheme, as described by
`Fixsen et al 2000 <https://ui.adsabs.harvard.edu/abs/2000PASP..112.1350F/abstract>`_.

Segments are determined using the 4-D GROUPDQ array of the input data set,
under the assumption that the :ref:`saturation detection <saturation_step>`
and :ref:`jump detection <jump_step>` steps have already been applied, in order
to flag occurrences of both saturation and cosmic-ray (CR) hits.
Ramps are broken into multiple segments where CR flags are found, so that slopes
are determined independently before and after the CR hits. Segments are
terminated at the first instance of a saturation flag.
Pixels are processed simultaneously in blocks using the array-based functionality of numpy.
The size of the block depends on the image size and the number of groups per
integration.

Upon successful completion of this step, the status keyword S_RAMP will be set
to "COMPLETE".

Output Products
---------------

RATE Product
++++++++++++
After computing the slopes and variances for all segments for a given pixel, the final slope is
determined as a weighted average from all segments in all integrations, and is
written to the "rate" output product.  In this output product, the
4-D GROUPDQ from all integrations is collapsed into 2-D, merged
(using a bitwise OR) with the input 2-D PIXELDQ, and stored as a 2-D DQ array. 
The 3-D VAR_POISSON and VAR_RNOISE arrays from all integrations are averaged
into corresponding 2-D output arrays.  In cases where the median rate
for a pixel is negative, the VAR_POISSON is set to zero, in order to avoid the
unphysical situation of having a negative variance.

RATEINTS Product
++++++++++++++++
The slope images for each integration are stored as a data cube in "rateints" output data
product.  Each plane of the 3-D SCI, ERR, DQ, VAR_POISSON, and VAR_RNOISE
arrays in this product corresponds to the result for a given integration.  In this output
product, the GROUPDQ data for a given integration is collapsed into 2-D and then
merged with the input 2-D PIXELDQ array to create the output DQ array for each
integration. The 3-D VAR_POISSON and VAR_RNOISE arrays are
calculated by averaging over the fit segments in the corresponding 4-D 
variance arrays.

FITOPT Product
++++++++++++++
A third, optional output product is also available and is produced only when
the step parameter `save_opt` is True (the default is False).  This optional
product contains 4-D arrays called SLOPE, SIGSLOPE, YINT, SIGYINT, WEIGHTS,
VAR_POISSON, and VAR_RNOISE, which contain the slopes, uncertainties in the
slopes, y-intercept, uncertainty in the y-intercept, fitting weights,
variance of the slope due to poisson noise, and the variance of the slope
due to read noise for each segment of each pixel, respectively. The y-intercept refers
to the result of the fit at an effective exposure time of zero.  This product also
contains a 3-D array called PEDESTAL, which gives the signal at zero exposure
time for each pixel, and the 4-D CRMAG array, which contains the magnitude of
each group that was flagged as having a CR hit.

By default, the name of this 
output file will have the product type suffix "_fitopt".
In this optional output product, the pedestal array is
calculated for each integration by extrapolating the final slope (the weighted
average of the slopes of all ramp segments in the integration) for each pixel
from its value at the first group to an exposure time of zero. Any pixel that is
saturated on the first group is given a pedestal value of 0. Before compression,
the cosmic-ray magnitude array is equivalent to the input SCI array but with the
only nonzero values being those whose pixel locations are flagged in the input
GROUPDQ as cosmic ray hits. The array is compressed, removing all groups in
which all the values are 0 for pixels having at least one group with a non-zero
magnitude. The order of the cosmic rays within the ramp is preserved.

Special Cases
-------------
If the input dataset has only one group in each integration (NGROUPS=1), the count rate
for all unsaturated pixels in each integration will be calculated as the
value of the science data in the one group divided by the group time.  If the
input dataset has only two groups per integration (NGROUPS=2), the count rate for all
unsaturated pixels in each integration will be calculated using the differences
between the two valid groups of the science data divided by the group time.

For datasets having more than one group in each integration (NGROUPS>1), a ramp having
a segment with only one good group is processed differently depending on the
number and size of the other segments in the ramp. If a ramp has only one
segment and that segment contains a single group, the count rate will be calculated
to be the value of the science data in that group divided by the group time.  If a ramp
has a segment with only one good group, and at least one other segment having more
than one good group, only data from the segment(s) having more than one
good group will be used to calculate the count rate.

For ramps in a given integration that are saturated beginning in their second group,
the count rate for that integration will be calculated as the value of the science data
in the first group divided by the group time, but only if the step parameter
``suppress_one_group`` is set to ``False``. If set to ``True``, the computation of
slopes for pixels that have only one good group will be suppressed and the slope
for that integration will be set to zero.

NIRCam Frame Zero
-----------------
The NIRCam instrument has the ability to downlink data resulting from the initial
frame of each integration (known as "frame zero") when on-board frame averaging is
in use for a given exposure. If the frame zero data were downlinked, they will appear
in a 3-D data cube in the raw data products with a label of "ZEROFRAME".
The data from frame zero can be used to recover a slope estimate for a pixel in the
event the pixel saturates already somewhere within the first *group* of an integration.

If all groups in an integration are flagged as SATURATED for a given pixel, the frame
zero data array is examined to determine whether or not it is also saturated. Saturated elements of
the frame zero array are set to zero by the preceding :ref:`saturation <saturation_step>`
step in the pipeline. Unsaturated elements will have non-zero values in the
frame zero array. If frame zero is *not* saturated, then it's value will be
divided by the frame time for the exposure in order to compute a slope for the pixel
in that integration. This is analagous to the situation in which only the first group
in an integration is unsaturated and used by itself to compute a slope.

Note that the computation of slopes from either a single group or the single frame
zero value is disabled when the step parameter ``suppress_one_group`` is set to ``True``.
In this case the slope value for such a pixel will be set to zero. <=== CHECK THIS.

Multiprocessing
---------------
This step has the option of running in multiprocessing mode. In that mode it will
split the input data into a number of slices based on the number of available
cores on the host computer and the value of the `maximum_cores` step parameter. By
default the step runs on a single processor. At the other extreme, if `maxiumum_cores` is
set to 'all', it will use all available cores (real and virtual). Testing has shown
a reduction in the elapsed time for the step proportional to the number of real
cores used. Using the virtual cores also reduces the elapsed time, but at a slightly
lower rate than the real cores.
Because the data are sliced based on the number
of rows, if the number of cores requested for multiprocessing is greater than
the number of rows, the number of cores actually used will be no more than the
number of rows.  This prevents any additional cores from operating on empty
datasets, which would cause errors during ramp fitting.

Detailed Algorithms
-------------------
The core algorithms for this step are called from the external package ``stcal``.
This package is an STScI effort to unify common calibration processing algorithms
for use by multiple observatories.

An in-depth discussion of the algorithms used to compute slopes and variances
is found :ref:`here <stcal:ramp_slopes_and_variances>`.

.. _ramp_error_propagation:

Error Propagation
-----------------
Error propagation in the ``ramp_fitting`` step is implemented by carrying along
the individual variances in the slope due to Poisson noise and read noise at all
levels of calculations. The total error estimate at each level is computed as
the square-root of the sum of the two variance estimates.

In each type of output product generated by the step, the variance in the slope
due to Poisson noise is stored in the "VAR_POISSON" extension, the variance in
the slope due to read noise is stored in the "VAR_RNOISE" extension, and the
total error is stored in the "ERR" extension. In the optional output product,
these arrays contain information for every segment used in the fitting for each
pixel. In the "rateints" product they contain values for each integration, and
in the "rate" product they contain values for the exposure as a whole.

.. _ramp_dq_propagation:

Data Quality Propagation
------------------------
For a given pixel, if all groups in an integration are flagged as DO_NOT_USE or
SATURATED, then that pixel will be flagged as DO_NOT_USE in the corresponding
integration in the "rateints" product.  Note this does NOT mean that all groups
are flagged as SATURATED, nor that all groups are flagged as DO_NOT_USE.  For
example, slope calculations that are suppressed due to a ramp containing only
one good group will be flagged as DO_NOT_USE in the
first group, but not necessarily any other group, while only groups two and
beyond are flagged as SATURATED.  Further, only if all integrations in the "rateints"
product are flagged as DO_NOT_USE, then the pixel will be flagged as DO_NOT_USE
in the "rate" product.

For a given pixel, if all groups in an integration are flagged as SATURATED,
then that pixel will be flagged as SATURATED and DO_NOT_USE in the corresponding
integration in the "rateints" product.  This is different from the above case in
that this is only for all groups flagged as SATURATED, not for some combination
of DO_NOT_USE and SATURATED.  Further, only if all integrations in the "rateints"
product are flagged as SATURATED, then the pixel will be flagged as SATURATED
and DO_NOT_USE in the "rate" product.

For a given pixel, if any group in an integration is flagged as JUMP_DET, then
that pixel will be flagged as JUMP_DET in the corresponding integration in the
"rateints" product.  That pixel will also be flagged as JUMP_DET in the "rate"
product.

.. _ramp_charge_migration:

Charge Migration Special Case
-----------------------------
If the :ref:`charge migration <charge_migration_step>`
step has been performed prior to ramp fitting, any group whose value exceeds the
``signal_threshold`` parameter value in that step will have been flagged with the
CHARGELOSS and DO_NOT_USE DQ flags. Due to the presence of the DO_NOT_USE flags,
such groups are excluded from all slope calculations.

It is still desired, however, to have a read noise variance value for such pixels
that is similar to pixels unaffected by charge migration, so an additional type of
variance is calculated, in which the groups flagged with CHARGELOSS are still included,
despite the fact that those groups do not get included in slope calculations.
This version of the readnoise variance is the one stored in the VAR_RNOISE extension
of the various output products from the step, so that it will be the one used later
in the pipeline flow in the :ref:`resample <resample_step>` step, if that step is
executed using Inverse Variance Map (IVM) weighting in the resampling process.

The original version of readnoise variance described earlier, where all groups flagged
with DO_NOT_USE are *not* included, is still used internally
in all other calculations involving readnoise variance.
