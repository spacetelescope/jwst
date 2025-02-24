.. _outlier-detection-imaging:

Imaging Data
============

This module serves as the interface for applying ``outlier_detection`` to direct
image observations, like those taken with MIRI, NIRCam, and NIRISS.
A :ref:`Stage 3 association <asn-level3-techspecs>`,
which is loaded into a :py:class:`~jwst.datamodels.ModelLibrary` object,
serves as the basic format for all processing performed by this step.
This routine performs the following operations:

#. Extract parameter settings for the input models and merge them with any user-provided values.
   See :ref:`outlier detection arguments <outlier_detection_step_args>` for the full list of parameters.

#. By default, resample all input images to the same output WCS. The resample process is
   controlled by the ``fillval``, ``pixfrac``, ``kernel``, and ``good_bits`` parameters;
   see the :ref:`outlier detection arguments <outlier_detection_step_args>` for more information.
   Resampling can be turned off with the ``resample_data`` parameter.

   * Compute an output WCS that is large enough to encompass all the input images.
   * Resample all images from the *same exposure* onto this output WCS to create a mosaic of all the detectors
     for that exposure.  This product is referred to as a "grouped mosaic" since it groups all the detectors
     from the same exposure into a single image. Each dither position will result in
     a separate grouped mosaic, so only a single exposure ever contributes to each pixel in these mosaics.
     An explanation of how all NIRCam multiple detector group mosaics are
     defined from `a single exposure or from a dithered set of exposures
     <https://jwst-docs.stsci.edu/near-infrared-camera/nircam-operations/nircam-dithers-and-mosaics>`_
     can be found here.
   * Fill in pixels that have no valid contribution from any input exposure with the value 
     specified by the ``fillval`` parameter.

#. If the ``save_intermediate_results`` parameter is set to True, write the resampled images to disk
   with the suffix ``_outlier_i2d.fits``. These are not saved if ``resample_data`` is set to False because
   they would be copies of the input models.

#. If resampling is turned off, use the input data itself in place of the resampled data
   for subsequent processing.

#. Construct and apply a bad pixel mask to the resampled data based on the drizzled weights.
   The ``weight_type`` parameter indicates the type of weighting image to apply with the bad pixel mask.
   The ``maskpt`` parameter sets the threshold weight value such that any pixel
   with a weight below this value gets flagged as bad.

#. Create a median image from all grouped observation mosaics pixel-by-pixel, i.e., averaging over groups.
   If ``save_intermediate_results`` is set to True, the median image is written out to disk with the
   suffix ``_median.fits``.

#. Blot (inverse of resampling) the median image back to match each original input image, and write 
   the blotted images to disk with the suffix ``_blot.fits`` if ``save_intermediate_results`` is `True`.

#. Perform statistical comparison between blotted image and original image to identify outliers.

   * If resampling is disabled (``resample_data == False``), compare the median image directly
     to each input image and ignore all sub-bullets below this one.
     In this case, compute the outlier mask using the following formula:

       .. math:: | image\_input - image\_median | > SNR * input\_err

   * Add a user-specified background value to the median image to match the original background levels
     of the input mosaic. This is controlled by the ``backg`` parameter.
   * Compute the spatial derivative of each pixel in the blotted median images by computing the absolute value
     of the difference between each pixel and its four surrounding neighbors, recording the largest
     absolute difference as the derivative. The derivative is multiplied by the ``scale`` parameter,
     allowing user flexibility in determining how aggressively to flag outliers at this stage.
   * Flag cosmic rays and other blemishes (such as satellite trails) based on the derivative image.
     Where the difference is larger than can be explained by a combination of noise statistics,
     the flattening effect of taking the median, and an error in the shift
     (the latter two effects are estimated using the image derivative), the suspect pixel is considered
     an outlier. The following rule is used:

     .. math:: | image\_input - image\_blotted | > scale*image\_deriv + SNR*noise

     The noise in this equation is the value of the input model's ``err`` extension.
     The SNR in this equation is the ``snr`` parameter, which encodes two values that
     determine whether a pixel should be masked:
     the first value detects the primary cosmic ray, and the second masks
     lower-level bad pixels adjacent to those found in the first pass. Since
     cosmic rays often extend across several pixels, the adjacent pixels make
     use of a slightly lower SNR threshold.

#. Update DQ arrays with flags and set SCI, ERR, and variance arrays to NaN at the location
   of identified outliers.

Memory Saving Options
---------------------
The outlier detection algorithm for imaging can require a lot of memory
depending on the number of inputs, the size of each input, and the size of the
final output product.  Specifically,

#. By default, all input exposures are kept open in memory to make
   processing more efficient.

#. The resample step creates an output product that is the
   same size as the final output product, which for imaging modes can span all detectors
   while also accounting for all dithers. Although only a single resampled image is needed in 
   memory at a time, for some Level 3 products, each resampled image can be on the order of several
   gigabytes in size.

#. The median combination step needs to have all pixels at the same position on
   the sky in memory in order to perform the median computation. The simplest (and fastest) implementation
   requires keeping all resampled outputs fully in memory at the same time.

These concerns have been addressed by implementing an overall memory model for outlier detection that
includes options to minimize memory usage at the expense of temporary file I/O and runtime.
Control over this memory model happens
with the use of the ``in_memory`` parameter, which defaults to True.
The full impact of setting this parameter to `False` includes:

#. The input :py:class:`~jwst.datamodels.ModelLibrary` object is loaded with `on_disk=True`.
   This ensures that input models are loaded into memory one at at time,
   and saved to a temporary file when not in use; these read-write operations are handled internally by
   the :py:class:`~jwst.datamodels.ModelLibrary` object.

#. Computing the median image works by writing the resampled data frames to appendable files
   on disk that are split into sections spatially but contain the entire ``groups``
   axis. The section size is set to use roughly the same amount of memory as a single resampled
   model, and since the resampled models are discarded from memory by the time the median calculation
   happens, this choice avoids increasing the overall memory usage of the step.
   Those sections are then read in one at a time to compute the median image.

These changes result in a minimum amount of memory usage during processing, but runtimes are
longer because many read and write operations are needed. Note that if a ModelLibrary object
is input to the step, the memory behavior of the step is read from the ``on_disk`` status
of the ModelLibrary object, and the ``in_memory`` parameter of the step is ignored.
When running ``calwebb_image3``, the ``in_memory`` flag should be set at the pipeline level,
e.g., ``strun calwebb_image3 asn.json --in-memory=False``; the step-specific flag will be ignored.
