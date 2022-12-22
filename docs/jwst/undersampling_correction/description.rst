Description
===========

:Class: `jwst.undersampling_correction.UndersamplingCorrectionStep`
:Alias: undersampling_correction

Overview
--------
This step corrects for an artifact seen in undersampled NIRISS images taken during commissioning.
The issue seen was that peak pixels of stars which reach close to saturation in dithers for which
the star is centered on the center of a pixel suffer from relatively strong charge migration such
that the group-to-group differences decrease significantly up the ramp once the signal level is
greater than roughly 30,000 ADU.  As a result, the last several groups of these ramps get flagged
by the ``jump`` step. The smaller number of groups used for these pixels in the ramp_fitting step
results in them having larger read noise variances, which in turn leads to lower weights used
during resampling. This ultimately leads to a lower than normal flux for the star in resampled
images. By flagging these groups affected by charge migration, they will not get used in the
computation of slopes in the ``ramp_fit`` step, and as detailed in the algorithm section below,
they will be used in the calculation of the variance of the slope due to readnoise,


Input details
-------------
The input data must have been processed through the ``jump`` step, so the input must be in the
form of a `~jwst.datamodels.RampModel`.


Algorithm
---------
The algorithm for this step is to flag as UNDERSAMP the groups exceeding the value of the
``signal_threshold`` parameter, which has units of ADU. Those groups will also be flagged as
DO_NOT_USE, so that they will not be included in the slope calculation in the ``ramp_fit``
step. The motivation for this step is to have the readnoise variance (calculated in ``ramp_fit``) be
similar to pixels unaffected by charge migration. Despite being also flagged as DO_NOT_USE, those
UNDERSAMP groups will be included in the calculation of the readnoise variance by ignoring the
DO_NOT_USE flag. For the Poisson noise variance calculation in ``ramp_fit``, UNDERSAMP/DO_NOT_USE
groups will not be included.

For integrations having only 1 or 2 groups, no flagging will be performed.


Output product
--------------
The output is a new copy of the input `~jwst.datamodels.RampModel`, with the corrections applied
to the GROUPDQ array.
