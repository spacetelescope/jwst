Description
===========

:Class: `jwst.emicorr.EmiCorrStep`
:Alias: emicorr

Overview
--------
The ``emicorr`` step corrects for known noise patterns in raw MIRI ramp data.
The majority of the MIRI subarrays have a 390 Hz or other electromagnetic
interference (EMI) noise pattern in the raw data.
The effect of these EMI frequencies is to imprint periodic noise into each ramp
image, depending on the readout timing of each pixel. For short integrations
in the LRS slitless mode, the correlated noise is quite apparent in the rate images.
For longer integrations the net effect is to increase the read noise by about 20\%.

In the long term, it may be possible to eliminate the most significant sources
of EMI noise by changing the sizes and locations of the subarrays, to get
readout times in phase with the known noise frequencies. For existing and
near-term observations, this step is required to correct for the noise.  In the
future, it is likely this step will still be useful to correct for lower level
EMI noise sources.

The known frequencies to correct for are stored in the
:ref:`EMICORR reference file <emicorr_reffile>`.
This reference file is expected to be in ASDF format, containing
the exact frequency numbers and the corresponding 1D array for the phase
amplitudes in the reference waveform.  The reference file also contains
a ``subarray_cases`` dictionary that maps a specific subarray, read pattern, and
detector to the necessary frequencies to correct for.

There are two algorithms available for fitting EMI noise in the input ramps: a
"sequential" fit and a "joint" fit.

For the sequential fitting algorithm, the process to correct for EMI noise is
the following (repeated iteratively for each discrete EMI frequency desired):

#. Read the image data.

#. Make a very crude slope image and fixed pattern bias image for each
   integration, ignoring the effects of nonlinearity, saturation, bad pixels, etc.

#. Subtract the scaled slope image and bias from each frame of each integration.

#. Calculate the phase of every pixel in the image at the desired EMI frequency
   (e.g. 390 Hz), relative to the first pixel in the image.

#. Make a binned phased amplitude (PA) waveform from the cleaned data.

#. Measure the phase shift between the binned PA waveform and the input
   reference waveform. [#f1]_

#. Use a look-up table to get the aligned reference waveform value for each pixel
   (the noise array corresponding to the input image). [#f1]_

#. Subtract the noise array from the input image and return the cleaned result.

.. [#f1] Alternately, if a reference waveform is not available, use the binned
   PA waveform to "self-correct" the noise.

The joint algorithm proceeds similarly, except that the linear ramps and
the EMI noise are fit simultaneously. It works by choosing pixels with modest
scatter among the reads, and then finding the amplitude and phase of a supplied
reference waveform that, when subtracted, makes these pixels' ramps as straight
as possible. The straightness of the ramps is measured by a chi-squared metric, after
fitting lines to each one.  As for the sequential algorithm, the EMI signal at each
frequency is fit and removed iteratively, for each desired frequency.

The optimal choice for the fitting algorithm depends on the input.
The sequential algorithm can be used without a reference waveform, generating
a new reference file on-the-fly for user-specified frequencies, but it
requires 10 or more groups for a reliable fit.  The joint algorithm
requires a reference waveform but can successfully fit EMI in ramps with
3 or more groups.

Input
-----
The input file is the ``_uncal`` file after the
:ref:`dq_init_step <dq_init_step>` step has been applied, in the
:ref:`calwebb_detector1 <calwebb_detector1>` pipeline.

Output
------
The output will be a partially-processed ``RampModel`` with the
EMI-corrected data in the SCI extension. All other extensions will
remain the same.
