.. _outlier_detection_:


Python Interface to OutlierDetection: OutlierDetection()
=========================================================

This module serves as the interface for applying outlier_detection to direct
image observations, like those taken with NIRCam.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST.

Specifically, this routine performs the following operations:

* Extract parameter settings from input model and merge them with any user-provided values

  - If resampling has been turned on for input CubeModel observations, convert the input CubeModel data into a ModelContainer.

<<<<<<< HEAD
* Resamples all input images into grouped observation mosaics; for example,
  combining all NIRCam multiple detector images from `a single exposure or
  from a dithered set of exposures.
  <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Dithers+and+Mosaics>`_
=======
* Resamples all input images into grouped observation mosaics.
>>>>>>> Update for outlier_detection documentation for B7.1

  - Resampled images will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

* Creates a median image from all grouped observation mosaics.

  - Median image will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

* Blot median image to match each original input exposure.

  - Resampled/blotted images will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

* Perform statistical comparison between blotted image and original image to identify outliers.
* Updates input data model DQ arrays with mask of detected outliers.


.. automodapi:: jwst.outlier_detection.outlier_detection
