.. _outlier-detection-imaging:

Default OutlierDetection Algorithm
-----------------------------------

This module serves as the interface for applying outlier_detection to direct
image observations, like those taken with NIRCam.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST.

Specifically, this routine performs the following operations:

* Extract parameter settings from input model and merge them with any user-provided values.
  The full set of user parameters includes::

    wht_type: type of data weighting to use during resampling;
              options are 'exptime','error','None' [default='exptime']
    pixfrac: pixel fraction used during resampling;
             valid values go from 0.0-1.0 [default=1.0]
    kernel: name of resampling kernel; options are 'square', 'turbo', 'point',
            'lanczos', 'tophat' [default='square']
    fillval: value to use to replace missing data when resampling;
              any floating point value (as a string) is valid (default='INDEF')
    nlow: Number (as an integer) of low values in each pixel stack to ignore
          when computing median value [default=0]
    nhigh:  Number (as an integer) of high values in each pixel stack to ignore
          when computing median value [default=0]
    maskpt: Percent of maximum weight to use as lower-limit for valid data;
            valid values go from 0.0-1.0 [default=0.7]
    grow: Radius (in pixels) from bad-pixel for neighbor rejection [default=1]
    snr: Signal-to-noise values to use for bad-pixel identification; valid
         values are a pair of floating-point values in a single string
         [default='4.0 3.0']
    scale: Scaling factor applied to derivative used to identify bad-pixels;
           valid value is a string with 2 floating point values [default='0.5 0.4')]
    backg: user-specified background value to apply to median image;
           [default=0.0]
    save_intermediate_results: specifies whether or not to save any products
                                created during outlier_detection [default=False]
    resample_data: specifies whether or not to resample the input data [default=True]
    good_bits: List of DQ integer values which should be considered good when
               creating weight and median images [default=0]

* Convert input data, as needed, to make sure it is in a format that can be processed

  - A :py:class:`~jwst.datamodels.ModelContainer` serves as the basic format for
    all processing performed by
    this step, as each entry will be treated as an element of a stack of images
    to be processed to identify bad-pixels/cosmic-rays and other artifacts.
  - If the input data is a :py:class:`~jwst.datamodels.CubeModel`, convert it into a ModelContainer.
    This allows each plane of the cube to be treated as a separate 2D image
    for resampling (if done at all) and for combining into a median image.

* By default, resample all input images into grouped observation mosaics; for example,
  combining all NIRCam multiple detector images from `a single exposure or
  from a dithered set of exposures.
  <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Dithers+and+Mosaics>`_

  - Resampled images will be written out to disk if
    ``save_intermediate_results`` parameter has been set to `True`
  - **If resampling was turned off**, a copy of the input (as a ModelContainer)
    will be used for subsequent processing.

* Create a median image from all grouped observation mosaics.

  - The median image will be created by combining all grouped mosaic images or
    non-resampled input data (as planes in a ModelContainer) pixel-by-pixel.
  - Median image will be written out to disk if ``save_intermediate_results``
    parameter has been set to `True`.

* By default, the median image will be blotted back (inverse of resampling) to
  match each original input exposure.

  - Resampled/blotted images will be written out to disk if
    ``save_intermediate_results`` parameter has been set to `True`
  - **If resampling was turned off**, the median image will be compared directly to
    each input image.
* Perform statistical comparison between blotted image and original image to identify outliers.
* Update input data model DQ arrays with mask of detected outliers.


.. automodapi:: jwst.outlier_detection.outlier_detection
