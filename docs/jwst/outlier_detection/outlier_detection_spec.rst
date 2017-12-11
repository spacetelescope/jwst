.. _outlier_detection_spec_:


Python Interface to OutlierDetectionSpec: OutlierDetectionSpec()
================================================================

This module serves as the interface for applying outlier_detection to long-slit
spectroscopic observations.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST
spectroscopic observations.

Specifically, this routine performs the following operations:

* Extract parameter settings from input model and merge them with any user-provided values

* Resamples all input images into a :py:class:`~jwst.datamodels.ModelContainer` using :py:class:`~jwst.resample.resample_spec.ResampleSpecData`

  - Resampled images will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

* Creates a median image from resampled :py:class:`~jwst.datamodels.ModelContainer`

  - Median image will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

* Blot median image to match each original input exposure.

  - Resampled/blotted images will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

* Perform statistical comparison between blotted image and original image to identify outliers.
* Updates input data model DQ arrays with mask of detected outliers.


.. automodapi:: jwst.outlier_detection.outlier_detection_spec
