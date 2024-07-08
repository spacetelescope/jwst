.. _outlier-detection-coron:

Coronography Outlier Detection Algorithm
========================================

This module serves as the interface for applying ``outlier_detection`` to coronographic
image observations.

Specifically, this routine performs the following operations:

#. Extract parameter settings from input model and merge them with any user-provided values.
   See :ref:`outlier detection arguments <outlier_detection_step_args>` for the full list
   of parameters.

#. Convert input data, as needed, to make sure it is in a format that can be processed.
   A :py:class:`~jwst.datamodels.CubeModel` serves as the basic format for all processing
   performed by this step.

#. Create a median image preserving the spatial dimensions of the cube.

   * The ``maskpt`` parameter sets the percentage of the weight image values to
     use, and any pixel with a weight below this value gets flagged as "bad" and
     ignored when resampled.

#. Perform statistical comparison between median image and original image to identify outliers.

   * This comparison uses the original input images, the
     median image, and the derivative of the median image to
     create a cosmic ray mask for each input image.
   * The derivative of the median image gets created using the
     median image to compute the absolute value of the difference between each pixel and
     its four surrounding neighbors with the largest value being the recorded derivative.
   * These derivative images are used to flag cosmic rays
     and other blemishes, such as satellite trails. Where the difference is larger
     than can be explained by noise statistics, the flattening effect of taking the
     median, or an error in the shift (the latter two effects are estimated using
     the image derivative), the suspect pixel is masked.
   * The ``backg`` parameter specifies a user-provided value to be used as the
     background estimate.  This gets added to the background-subtracted
     median image to attempt to match the original background levels of the
     original input mosaic so that cosmic-rays (bad pixels) from the input
     mosaic can be identified more easily as outliers compared to the median image.
   * Cosmic rays are flagged using the following rule:

     .. math:: | image\_input - image\_median | > scale*image\_deriv + SNR*noise

   * The ``scale`` is defined as the multiplicative factor applied to the
     derivative which is used to determine if the difference between the data
     image and the median image is large enough to require masking.
   * The ``noise`` is calculated using a combination of the detector read
     noise and the poisson noise of the median image plus the sky background.
   * The user must specify two cut-off signal-to-noise values using the
     ``snr`` parameter for determining whether a pixel should be masked:
     the first for detecting the primary cosmic ray, and the second for masking
     lower-level bad pixels adjacent to those found in the first pass. Since
     cosmic rays often extend across several pixels, the adjacent pixels make
     use of a slightly lower SNR threshold.

#. Update input data model DQ arrays with mask of detected outliers.
