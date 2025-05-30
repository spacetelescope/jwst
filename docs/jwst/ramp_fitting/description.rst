Description
===========

:Class: `jwst.ramp_fitting.RampFitStep`
:Alias: ramp_fit

This step determines the mean count rate, in units of counts per second, for
each pixel by performing a linear fit to the "up the ramp" data in the input file.
All fitting is done using the "ordinary least squares" (OLS_C) method.
The fit is performed independently for each pixel. Bad values flagged via
DQ flags are rejected from the fits.

The input to the step is assumed to be the fully-corrected and flagged 4-D
data resulting from applying all previous steps of the
:ref:`calwebb_detector1 <calwebb_detector1>` pipeline and will nominally be
the output from the :ref:`jump detection <jump_step>` step. It is in fact
vital that all anomalies such as saturation, non-linearity, and CR jumps
be corrected or appropriately flagged in order to obtain useful results
from ramp fitting.

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
A ramp segment is a set of contiguous groups that have no non-zero DQ values
assigned. The one exception to this rule is the occurrence of a "JUMP_DET"
(jump detected) flag: a group with this flag will be used as the first group of
the next segment. Any occurrences of a "DO_NOT_USE" flag will be excluded from a
segment. When a "SATURATION" flag is found, the segment is terminated at the
preceding group and all subsequent groups are rejected.
Any segment containing only one good group is ignored if there is any other
segment of length greater than one.
Once all segments have been determined, slopes and variances are determined for
each one.

The range of groups to be fitted can be limited using the ``--firstgroup`` and
``--lastgroup`` parameters.  This works by setting the DO_NOT_USE DQ bit in the GROUPDQ
attribute for all groups outside the range selected.

Pixels are processed simultaneously in blocks using the array-based functionality of numpy.
The size of the block depends on the image size and the number of groups per
integration.

The main algorithms for this step are called from the external package ``stcal``.
This package is an STScI effort to unify common calibration processing algorithms
for use by multiple observatories.
Therefore the majority of the remainder of this document links to the relevant
sections of information in the ``stcal`` package.
JWST-specific features are described later within this document.

Upon successful completion of this step, the status keyword S_RAMP will be set
to "COMPLETE".

There is a new algorithm available for testing in ramp fitting, which is the
likelihood algorithm.  It is selected by setting ``--ramp_fitting.algorithm=LIKELY``.
The details are in the ``stcal`` documentation at
:ref:`Likelihood Algorithm Details <stcal:likelihood_algo>`.

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

:ref:`Output Products <stcal:ramp_output_products>`
---------------------------------------------------

:ref:`Special Cases <stcal:ramp_special_cases>`
-----------------------------------------------

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
in that integration. This is analogous to the situation in which only the first group
in an integration is unsaturated and used by itself to compute a slope.

Note that the computation of slopes from either a single group or the single frame
zero value is disabled when the step parameter ``suppress_one_group`` is set to ``True``.
In this case the slope value for such a pixel will be set to zero.

:ref:`Detailed Algorithms <stcal:ramp_slopes_and_variances>`
------------------------------------------------------------

:ref:`Error Propagation <stcal:ramp_error_propagation>`
-------------------------------------------------------

:ref:`Data Quality Propagation <stcal:ramp_dq_propagation>`
-----------------------------------------------------------

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
