Description
===========

:Class: `jwst.undersampling_correction.UndersamplingCorrectionStep'
:Alias: undersampling_correction

Overview
--------
This step corrects for an artifact seen in undersampled NIRISS images taken during commissioning.
The issue seen was that peak pixels of stars which reach close to saturation in dithers for which
the star is centered on the center of a pixel suffer from relatively strong charge migration such
that the group-to-group differences decrease significantly up the ramp once the signal level is
greater than roughly 30,000 ADU.  As a result, the last several groups of these ramps get flagged
by the ``jump`` step, even though those jumps are negative.  Consequently, those groups are
assigned low weights during the ``ramp_fitting`` and ``resampling`` steps, resulting in a
relatively low integrated flux for the star.


Input details
-------------
The input data must have been processed through the ``jump`` step, so the input must be in the
form of a `~jwst.datamodels.RampModel`.


Algorithm
---------
The algorithm for this step is to flag as DO_NOT_USE the groups exceeding the value of the
``signal_threshold`` parameter, which has units of ADU.


Output product
--------------
The output is a new copy of the input `~jwst.datamodels.RampModel`, with the corrections applied
to the GROUPDQ array.
