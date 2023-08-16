Description
===========

:Class: `jwst.undersampling_correction.UndersamplingCorrectionStep`
:Alias: undersampling_correction

Overview
--------
This step corrects for an artifact seen in undersampled NIRISS images that may depress flux 
in resampled images. The artifact is seen in dithered images where the star is centered in 
a pixel. When the peak pixels of such stars approach the saturation level, they suffer from 
significant :ref:`charge migration <charge_migration>`:
the spilling of charge from a saturated pixel into its neighboring pixels. This charge migration 
causes group-to-group differences to decrease significantly once the signal level is greater than 
~30,000 ADU. As a result, the last several groups of these ramps get flagged by the ``jump`` step. 
The smaller number of groups used for these pixels in the ``ramp_fitting`` step results in them having 
larger read noise variances, which in turn leads to lower weights used during resampling. This 
ultimately leads to a lower than normal flux for the star in resampled images.

Once a group in a ramp has been flagged as affected by charge migration, all subsequent 
groups in the ramp are also flagged. By flagging these groups, they are not used in the 
computation of slopes in the :ref:`ramp_fitting <ramp_fitting_step>` step. However, as described 
in the algorithm section below, they _are_ used in the calculation of the variance of the slope 
due to readnoise.

Input details
-------------
The input data must have been processed through the ``jump`` step, so the input must be in the
form of a `~jwst.datamodels.RampModel`.


Algorithm
--------- 
The first group, and all subsequent groups, exceeding the value of the 
``signal_threshold`` parameter is flagged as UNDERSAMP. ``signal_threshold`` is in units 
of ADUs. These groups will also be flagged as DO_NOT_USE, and will not 
be included in the slope calculation during the ``ramp_fitting`` step. Despite being flagged 
as DO_NOT_USE, these UNDERSAMP groups are still included in the calculation of the
variance due to readnoise. 
This results in a readnoise variance for undersampled pixels that is similar to that of 
pixels unaffected by charge migration. For the Poisson noise variance calculation in 
:ref:`ramp_fitting <ramp_fitting_step>`, the UNDERSAMP/DO_NOT_USE groups are not included.

For integrations having only 1 or 2 groups, no flagging will be performed.


Output product
--------------
The output is a new copy of the input `~jwst.datamodels.RampModel`, with the updated DQ flags
added to the GROUPDQ array.

