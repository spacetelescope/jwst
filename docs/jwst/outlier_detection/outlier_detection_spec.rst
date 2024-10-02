.. _outlier-detection-spec:

Outlier Detection for Slit-like Spectroscopic Data
==================================================

This module serves as the interface for applying ``outlier_detection`` to slit-like
spectroscopic observations.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST
spectroscopic observations.

Specifically, this routine performs the following operations (modified from the
:ref:`Default Outlier Detection Algorithm <outlier-detection-imaging>` ):

#. Extract parameter settings from input model and merge them with any user-provided values

   - The same set of parameters available to:
     ref:`Default Outlier Detection Algorithm <outlier-detection-imaging>`
     also applies to this code.

#. Convert input data, as needed, to make sure it is in a format that can be processed

   - A :py:class:`~jwst.datamodels.ModelContainer` serves as the basic format 
     for all processing performed by
     this step, as each entry will be treated as an element of a stack of images
     to be processed to identify bad pixels, cosmic-rays and other artifacts.

#. If the ``resample_data`` parameter is set to True, resample all input images
   using :py:class:`~jwst.resample.resample_spec.ResampleSpecData`.

   - Error images are resampled alongside the science data, to create
     approximate error arrays for each resampled exposure.
   - The resampled data are written out to disk with suffix "outlier_s2d"
     if the ``save_intermediate_results`` parameter is set to `True`.

#. Create a median image from the resampled exposures, or directly from
   the input exposures, if ``resample_data`` is set to False.

   - The error images for each exposure are also median-combined.
   - The median data are written out to disk with suffix "median"
     if the ``save_intermediate_results`` parameter is set to `True`.

#. If ``resample_data`` is set to True, blot the median image to match
   each original input image.

   - The median error image is also blotted to match the input.
   - Resampled/blotted images are written out to disk if the ``save_intermediate_results``
     parameter is set to `True`.

#. Perform statistical comparison between the median image and the original image
   to identify outliers.

   - Signal-to-noise thresholding uses the median error image, rather than the
     input error image, in order to better identify outliers that have both
     high data values and high error values in the input exposures.

#. Update input data model DQ arrays with mask of detected outliers.
   Update the SCI, ERR, and VAR arrays with NaN values at the outlier
   locations.


.. automodapi:: jwst.outlier_detection.spec
