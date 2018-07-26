.. _outlier-detection-spec:

OutlierDetection for Long-Slit Spectroscopic Data
====================================================

This module serves as the interface for applying outlier_detection to long-slit
spectroscopic observations.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST
spectroscopic observations.

Specifically, this routine performs the following operations (modified from the
:ref:`Default Outlier Detection Algorithm <outlier-detection-imaging>` ):

* Extract parameter settings from input model and merge them with any user-provided values

  - the same set of parameters available to :ref:`Default Outlier Detection Algorithm <outlier-detection-imaging>`
    also applies to this code

* Convert input data, as needed, to make sure it is in a format that can be processed

  - A :py:class:`~jwst.datamodels.ModelContainer` serves as the basic format 
    for all processing performed by
    this step, as each entry will be treated as an element of a stack of images
    to be processed to identify bad-pixels/cosmic-rays and other artifacts.
  - If the input data is a :py:class:`~jwst.datamodels.CubeModel`, convert it into a 
    :py:class:`~jwst.datamodels.ModelContainer`.
    This allows each plane of the cube to be treated as a separate 2D image
    for resampling (if done at all) and for combining into a median image.

* Resamples all input images into a :py:class:`~jwst.datamodels.ModelContainer` using :py:class:`~jwst.resample.resample_spec.ResampleSpecData`

  - Resampled images will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`
  - **If resampling was turned off**, the original inputs will be used to create
    the median image for cosmic-ray detection.

* Creates a median image from (possibly) resampled :py:class:`~jwst.datamodels.ModelContainer`

  - Median image will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`

* Blot median image to match each original input exposure.

  - Resampled/blotted images will be written out to disk if ``save_intermediate_results`` parameter has been set to `True`
  - **If resampling was turned off**, the median image will be used as for comparison
    with the original input models for detecting cosmic-rays.

* Perform statistical comparison between blotted image and original image to identify outliers.
* Updates input data model DQ arrays with mask of detected outliers.


.. automodapi:: jwst.outlier_detection.outlier_detection_spec
