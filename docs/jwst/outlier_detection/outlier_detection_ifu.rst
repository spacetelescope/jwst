.. _outlier_detection_ifu_:


Python Interface to OutlierDetectionIFU: OutlierDetectionIFU()
=========================================================

This module serves as the interface for applying outlier_detection to IFU
 observations, like those taken with NIRSpec and MIRI.  The code implements the
 basic outlier detection algorithm used with HST data, as adapted to JWST IFU
 observations.

 Specifically, this routine performs the following operations:

 * Extract parameter settings from input model and merge them with any user-provided values

 * Resamples all input images into :py:class:`~jwst.datamodels.IFUCubeModel` observations.

   - Resampled cubes will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

 * Creates a median image from the set of resampled :py:class:`~jwst.datamodels.IFUCubeModel` observations

   - Median image will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

 * Blot median image to match each original input exposure.

   - Resampled/blotted cubes will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

 * Perform statistical comparison between blotted image and original image to identify outliers.
 * Updates input data model DQ arrays with mask of detected outliers.


.. automodapi:: jwst.outlier_detection.outlier_detection_ifu
