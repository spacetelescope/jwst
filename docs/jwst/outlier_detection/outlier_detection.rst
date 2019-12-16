.. _outlier-detection-imaging:

Default OutlierDetection Algorithm
===================================

This module serves as the interface for applying ``outlier_detection`` to direct
image observations, like those taken with NIRCam and MIRI.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST.

Specifically, this routine performs the following operations:

* Extract parameter settings from input model and merge them with any user-provided values.
  The full set of user parameters includes::

    weight_type: type of data weighting to use during resampling;
                 options are 'exptime', 'error', 'None' [default='exptime']
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
    good_bits: Sum of DQ integer values which should be considered good when
               creating weight and median images [default=6]
    scale_detection: Boolean indicating whether to rescale the individual input
                     images/integrations to match total signal when doing
                     comparisons [default=False]

* Convert input data, as needed, to make sure it is in a format that can be processed

  - A :py:class:`~jwst.datamodels.ModelContainer` serves as the basic format for
    all processing performed by
    this step, as each entry will be treated as an element of a stack of images
    to be processed to identify bad-pixels/cosmic-rays and other artifacts.
  - If the input data is a :py:class:`~jwst.datamodels.CubeModel`, convert it into a ModelContainer.
    This allows each plane of the cube to be treated as a separate 2D image
    for resampling (if done) and for combining into a median image.

* By default, resample all input images into grouped observation mosaics; for example,
  combining all NIRCam multiple detector images from `a single exposure or
  from a dithered set of exposures.
  <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Dithers+and+Mosaics>`_

  - Resampled images will be written out to disk if the
    ``save_intermediate_results`` parameter is set to `True`
  - **If resampling is turned off**, a copy of the input (as a ModelContainer)
    will be used for subsequent processing.

* Create a median image from all grouped observation mosaics.

  - The median image is created by combining all grouped mosaic images or
    non-resampled input data (as planes in a ModelContainer) pixel-by-pixel.
  - The median image is written out to disk if the ``save_intermediate_results``
    parameter is set to `True`.

* By default, the median image is blotted back (inverse of resampling) to
  match each original input image.

  - Resampled/blotted images are written out to disk if
    the ``save_intermediate_results`` parameter is set to `True`
  - **If resampling is turned off**, the median image is compared directly to
    each input image.
* Perform statistical comparison between blotted image and original image to identify outliers.
* Update input data model DQ arrays with mask of detected outliers.

Outlier Detection for TSO data
-------------------------------
Time-series observations (TSO) result in input data stored as a 3D CubeModel
where each plane in the cube represents a separate integration without changing the
pointing.  Normal imaging data would benefit from combining all integrations into a
single image. TSO data's value, however, comes from looking for variations from one
integration to the next.  The outlier detection algorithm, therefore, gets run with 
a few variations to accomodate the nature of these 3D data.

* Input data is converted from a CubeModel (3D data array) to a ModelContainer

  - Each model in the ModelContainer is a separate plane from the input CubeModel

* The median image is created without resampling the input data

  - All integrations are aligned already, so no resampling or shifting needs to be performed
  
* A matched median gets created by combining the single median frame with the 
  noise model for each input integration.

* Perform statistical comparison between the matched median with each input 
  integration.  

* Update input data model DQ arrays with the mask of detected outliers


.. note:: 

  This same set of steps also gets used to perform outlier detection on
  coronographic data, because it too is processed as 3D (per-integration)
  cubes.


Outlier Detection for IFU data
------------------------------
Integral Field Unit (IFU) data is handled as a 2D image on input (i.e. the state of
the data before creating a 3D cube).  This 2D image
gets converted into a properly calibrated spectral cube (3D array) and stored as
an IFUCubeModel for use within outlier detection.  The many differences in data format 
for the IFU data relative to normal direct imaging data requires special 
processing in order to perform outlier detection on IFU data.  

* Convert the input 2D IFUImageModel into a 3D IFUCubeModel by calling
  :py:class:`~jwst.cube_build.CubeBuildStep`

  - A separate IFUCubeModel is generated for each wavelength channel/band by
    using  the `single` option for the :py:class:`~jwst.cube_build.CubeBuildStep`.
    
* All IFUCubeModels get median combined to create a single median 
  IFUCubeModel product.
  
* The IFUCubeModel median product gets resampled back to match each 
  original input IFUImageModel dataset.
  
  - This resampling uses :py:class:`~jwst.cube_build.blot_cube_build.CubeBlot` 
    to perform this conversion.

* The blotted, median 2D images are compared statistically to the original 
  2D input images to detect outliers.
  
* The DQ array of each input dataset gets updated to mark the detected
  outliers.
  

.. automodapi:: jwst.outlier_detection.outlier_detection
