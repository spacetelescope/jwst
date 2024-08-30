.. _outlier-detection-coron:

Outlier Detection for Coronagraphic Data
========================================

This module serves as the interface for applying ``outlier_detection`` to coronagraphic
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

   The core detection algorithm uses the following to generate an outlier mask

   .. math:: | image\_input - image\_median | > SNR*input\_err

#. Update input data model DQ arrays with mask of detected outliers.

.. automodapi:: jwst.outlier_detection.coron
