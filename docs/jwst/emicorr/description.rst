Description
===========

:Class: `jwst.emicorr.EmiCorrStep`
:Alias: emicorr

Overview
--------
The ``emicorr`` step corrects for known noise patterns in the raw MIRI data.
The majority of the MIRI subarrays have an 390 Hz or other electromagnetic
interference (EMI) noise pattern in the raw data. The known frequencies to
correct for are in the EMI reference file, under the key ``frequencies``.
The effect of these EMI frequencies is to imprint each into the rate
images. For short integrations in ``LRSSLITLESS`` the correlated noise from
this is quite apparent in the rate images. For longer integrations the net
effect is to increase the read noise by about 20\%.

The process to do the correction is the following (repeated
recursively for each discrete EMI frequency desired):

#. Read image data.

#. Make very crude slope image and fixed pattern "super" bias for each
   integration, ignoring everything (nonlin, saturation, badpix, etc).

#. Subtract scaled slope image and bias from each frame of each integration.

#. Calculate phase of every pixel in the image at the desired EMI frequency
   (e.g. 390 Hz) relative to the first pixel in the image.

#. Make a binned, phased amplitude (pa) wave from the cleaned data (plot
   cleaned vs phase, then bin by phase).

#. Measure the phase shift between the binned pa wave and the input
   reference wave. [#f1]_

#. Use look-up table to get the aligned reference wave value for each pixel
   (the noise array corresponding to the input image). [#f1]_

#. Subtract the noise array from the input image and return the cleaned result.

.. [#f1] Alternately, use the binned pa wave instead of the reference wave to "self-correct".

The long term plan is a change to the sizes and locations of the subarrays
to get the frame times to be in phase with the known noise frequencies like
the full frame images. For the previous and near term observations this can
be fixed through application of the ``emicorr`` step.

An EMICORR reference file can be used to correct for all known noise
patterns. The reference file is expected to be in ASDF format, containing
the exact frequency numbers, the corresponding 1D array for the phase
amplitudes, and a ``subarray_cases`` dictionary that contains
the frequencies to correct for according to subarray, read pattern, and
detector. If there is no reference file found in CRDS, the step has a set
of default frequencies and subarray cases for which the correction is
applied.

Input
-----
The input file is the ``_uncal`` file after the
:ref:`dq_init_step <dq_init_step>` step has been
applied, in the in the :ref:`calwebb_detector1 <calwebb_detector1>`
pipeline.

Output
------
The output will be a partially-processed ``RampModel`` with the
corrected data in the SCI extension, meaning, the effect of the
EMI frequencies (either the default values or the ones in the
reference file) removed. All other extensions will remain the same.
