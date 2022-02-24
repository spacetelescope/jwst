.. _outlier-detection-ifu:

Outlier Detection for IFU Data
==============================

This module serves as the interface for applying ``outlier_detection`` to IFU
observations, like those taken with NIRSpec and MIRI.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST IFU
observations.

Specifically, this routine performs the following operations (modified from 
:ref:`Default Outlier Detection Algorithm <outlier-detection-imaging>` ):

 * Extract parameter settings from input model and merge them with any user-provided values

   - the same set of parameters available to 
     :ref:`Default Outlier Detection Algorithm <outlier-detection-imaging>` 
     also applies to this code

 * Resample all input :py:class:`~jwst.datamodels.IFUImageModel` images into 
   :py:class:`~jwst.datamodels.IFUCubeModel` data cubes

   - Resampling uses :py:class:`~jwst.cube_build.CubeBuildStep` to create 
     :py:class:`~jwst.datamodels.IFUCubeModel` formatted data for processing
   - Resampled cubes are written out to disk if the ``save_intermediate_results`` 
     parameter is set to `True`

 * Create a median cube from the set of resampled 
   :py:class:`~jwst.datamodels.IFUCubeModel` cubes

   - The median cube is written out to disk if the ``save_intermediate_results`` 
     parameter is set to `True`

 * Blot median cube to match each original 2D input image

   - Resampled/blotted images are written out to disk if 
     the ``save_intermediate_results`` parameter is set to `True`

 * Perform statistical comparison between blotted image and original image to identify outliers
 * Update input data model DQ arrays with mask of detected outliers


.. automodapi:: jwst.outlier_detection.outlier_detection_ifu
