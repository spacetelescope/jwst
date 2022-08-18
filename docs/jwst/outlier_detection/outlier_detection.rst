.. _outlier-detection-imaging:

Default Outlier Detection Algorithm
===================================

This module serves as the interface for applying ``outlier_detection`` to direct
image observations, like those taken with NIRCam and MIRI.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST.

Specifically, this routine performs the following operations:

* Extract parameter settings from input model and merge them with any user-provided values.
  See :ref:`outlier detection arguments <outlier_detection_step_args>` for the full list
  of parameters.

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
  <https://jwst-docs.stsci.edu/near-infrared-camera/nircam-operations/nircam-dithers-and-mosaics>`_

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

Memory Model for Outlier Detection Algorithm
---------------------------------------------
The outlier detection algorithm can end up using massive amounts of memory
depending on the number of inputs, the size of each input, and the size of the
final output product.  Specifically,

    * The input :py:class:`~jwst.datamodels.ModelContainer` or
      :py:class:`~jwst.datamodels.CubeModel`
      for IFU data, by default, all input exposures would have been kept open in memory to make
      processing more efficient.

    * The initial resample step creates an output product for EACH input that is the
      same size as the final
      output product, which for imaging modes can span all chips in the detector while
      also accounting for all dithers.  For some Level 3 products, each resampled image can
      be on the order of 2Gb or more.

    * The median combination step then needs to have all pixels at the same position on
      the sky in memory in order to perform the median computation.  The simplest implementation
      for this step requires keeping all resampled outputs fully in memory at the same time.

Many Level 3 products only include a modest number of input exposures which can be
processed using less than 32Gb of memory at a time.  However, there are a number of
ways this memory limit can be exceeded.  This has been addressed by implementing an
overall memory model for the outlier detection that includes options to minimize the
memory usage at the expense of file I/O.  The control over this memory model happens
with the use of the ``in_memory`` parameter.  The full impact of this parameter
during processing includes:

    * The ``save_open`` parameter gets set to `False`
      when opening the input :py:class:`~jwst.datamodels.ModelContainer` object.
      This forces all input models in the input :py:class:`~jwst.datamodels.ModelContainer` or
      :py:class:`~jwst.datamodels.CubeModel` to get written out to disk.  The ModelContainer
      then uses the filename of the input model during subsequent processing.

    * The ``in_memory`` parameter gets passed to the :py:class:`~jwst.resample.ResampleStep`
      to set whether or not to keep the resampled images in memory or not.  By default,
      the outlier detection processing sets this parameter to `False` so that each resampled
      image gets written out to disk.

    * Computing the median image works section-by-section by only keeping 1Mb of each input
      in memory at a time.  As a result, only the final output product array for the final
      median image along with a stack of 1Mb image sections are kept in memory.

    * The final resampling step also avoids keeping all inputs in memory by only reading
      each input into memory 1 at a time as it gets resampled onto the final output product.

These changes result in a minimum amount of memory usage during processing at the obvious
expense of reading and writing the products from disk.


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
