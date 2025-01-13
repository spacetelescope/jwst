.. _outlier-detection-coron:

Overview
========

This module serves as the interface for detecting outliers in coronagraphic
image observations, namely those taken in the following modes:
'MIR_LYOT', 'MIR_4QPM', 'NRC_CORON'.

This routine performs the following operations:

#. Extract parameter settings from input model and merge them with any user-provided values.
   See :ref:`arguments <outlier_detection_coron_step_args>` for the full list
   of parameters.

#. Do not attempt resampling; data are assumed to be aligned and have an identical WCS.
   This is true automatically for a CubeModel.

#. Create a median image over the `groups` (exposures, planes of cube) axis,
   preserving the spatial (x,y) dimensions of the cube.

   * The ``maskpt`` parameter sets the percentage of the weight image values to
     use, and any pixel with a weight below this value gets flagged as "bad".

#. Perform statistical comparison between median image and original image to identify outliers.

   The core detection algorithm uses the following to generate an outlier mask

   .. math:: | image\_input - image\_median | > SNR*input\_err

#. Update DQ arrays with flags and set SCI, ERR, and variance arrays to NaN at the location
   of identified outliers.

``OutlierDetectionCoronStep`` does not use any reference files.
