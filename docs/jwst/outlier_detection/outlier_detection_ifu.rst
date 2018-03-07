.. _outlier_detection_ifu_:


Python Interface to OutlierDetectionIFU: OutlierDetectionIFU()
=========================================================

This module serves as the interface for applying outlier_detection to IFU
<<<<<<< HEAD
 observations, like those taken with NIRSpec and MIRI.  The code implements the
=======
 observations, like those taken with NIRISS and MIRI.  The code implements the
>>>>>>> Update for outlier_detection documentation for B7.1
 basic outlier detection algorithm used with HST data, as adapted to JWST IFU
 observations.

 Specifically, this routine performs the following operations:

 * Extract parameter settings from input model and merge them with any user-provided values

 * Resamples all input images into :py:class:`~jwst.datamodels.IFUCubeModel` observations.

<<<<<<< HEAD
   - Resampled cubes will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`
=======
   - Resampled images will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`
>>>>>>> Update for outlier_detection documentation for B7.1

 * Creates a median image from resampled :py:class:`~jwst.datamodels.IFUCubeModel` observations

   - Median image will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

 * Blot median image to match each original input exposure.

<<<<<<< HEAD
   - Resampled/blotted cubes will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`
=======
   - Resampled/blotted images will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`
>>>>>>> Update for outlier_detection documentation for B7.1

 * Perform statistical comparison between blotted image and original image to identify outliers.
 * Updates input data model DQ arrays with mask of detected outliers.


.. automodapi:: jwst.outlier_detection.outlier_detection_ifu
