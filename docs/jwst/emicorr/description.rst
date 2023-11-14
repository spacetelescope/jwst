Description
===========

:Class: `jwst.emicorr.EmiCorrStep`
:Alias: emicorr

Overview
--------
The ``emicorr`` step corrects for known noise patterns in the raw MIRI data.
The majority of the MIRI subarrays have an 390 Hz or other electromagnetic
interference (EMI) noise pattern in the raw data. The known frequencies to
correct for are \[\Hz\]\: 390.625, 218.52055, 218.520470, 218.520665, 164.9305,
10.039216. For 390 Hz, the wave maps exactly to 256 pixels and has a constant
amplitude of \+\- 4 DN. The effect of this is to imprint this into the rate
images. For short integrations in LRSSLITLESS the correlated noise from this
is quite apparent in the rate images. For longer integrations the net effect
is to increase the read noise by about 20\%\.

The process to do the correction is the following (repeated
recursively for each discrete EMI frequency desired):
1. Read image data.
2. Make very crude slope image and fixed pattern "super" bias for each
integration, ignoring everything (nonlin, saturation, badpix, etc).
3. Subtract scaled slope image and bias from each frame of each integration.
4. Calculate phase of every pixel in the image at the desired EMI frequency
(e.g. 390 Hz) relative to the first pixel in the image.
5. Make a binned, phased amplitude (pa) wave from the cleaned data (plot
cleaned vs phase, then bin by phase).
6.* Measure the phase shift between the binned pa wave and the input
reference wave
7.* Use look-up table to get the aligned reference wave value for each pixel
(the noise array corresponding to the input image).
*) Alternately, use the binned pa wave instead of the reference wave to
"self-correct"
8. Subtract the noise array from the input image and return the cleaned result.

The long term plan is a change to the sizes and locations of the subarrays
to get the frame times to be in phase with the known noise frequencies like
the full frame images. For the previous and near term observations this can
be fixed the ``emicorr`` step.

An EMICORR reference file can be used for to correct for all known noise
patterns. The reference file is expected to be in an ASDF format, containing
the phase amplitude values corresponding to the known frequencies.

Input
-----
The input file is the ``_uncal`` file after the ``dq_init`` step has been
ran, in the in the calwebb_detector1 pipeline.

Output
------
The output will be an ``_uncal`` file with the corrected data in the SCI
extension, meaning, the effect of the known EMI frequencies subtracted. All
other extensions will remain the same.
