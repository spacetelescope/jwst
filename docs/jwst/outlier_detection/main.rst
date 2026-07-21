.. _outlier_design:

Overview
========

:Class: `jwst.outlier_detection.outlier_detection_step.OutlierDetectionStep`
:Alias: outlier_detection

This step provides the sole interface to all methods of performing outlier
detection on JWST observations.

Processing multiple datasets together allows for the identification of bad pixels
or cosmic rays that remain in each of the input images, often at levels which
were not detectable by the :ref:`jump <jump_step>` step.
The ``outlier_detection`` step supports multiple
algorithms and determines the appropriate algorithm for the type of observation
being processed.  This step supports:

* :ref:`Imaging modes <outlier-detection-imaging>`
    - Exposure types: 'FGS_IMAGE', 'MIR_IMAGE', 'NRC_IMAGE', 'NIS_IMAGE'
* :ref:`Slit-like Spectroscopic modes <outlier-detection-spec>`
    - Exposure types: 'MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT', 'NRS_MSASPEC'
* :ref:`Time-Series-Observation (TSO) modes <outlier-detection-tso>`
    - Exposure types: 'MIR_LRS-SLITLESS', 'NRC_TSGRISM', 'NIS_SOSS', 'NRS_BRIGHTOBJ',
      and 'NRC_TSIMAGE', as well as TSOs obtained with the MIRI imager or in fixed-slit mode
      ('MIR_IMAGE' or 'MIR_LRS-FIXEDSLIT' with ``TSOVISIT=True``)
* :ref:`IFU Spectroscopic modes <outlier-detection-ifu>`
    - Exposure types: 'MIR_MRS', 'NRS_IFU'
* :ref:`Coronagraphic Image modes <outlier-detection-coron>`
    - Exposure types: 'MIR_LYOT', 'MIR_4QPM', 'NRC_CORON'

This step uses the following logic to apply the appropriate algorithm to the
input data:

#. Interpret inputs (`~jwst.associations.Association`,
   `~jwst.datamodels.container.ModelContainer`,
   `~jwst.datamodels.library.ModelLibrary`, or
   `~stdatamodels.jwst.datamodels.CubeModel`)
   to identify all input observations to be processed.

#. Read in parameters set by user. See :ref:`outlier_detection_step_args` for the full list
   of parameters.

#. Select outlier detection algorithm based on exposure type in input model ``meta.exposure.type``.

#. Instantiate and run outlier detection class determined for the exposure type
   using parameter values interpreted from inputs.

#. Update DQ arrays with flags and set SCI, ERR, and variance arrays to NaN at the location
   of identified outliers.

.. _outlier-detection-imaging:

Imaging Data
------------

This section explains the interface for applying ``outlier_detection`` to direct
image observations, like those taken with MIRI, NIRCam, and NIRISS.
A :ref:`Stage 3 association <asn-level3-techspecs>`,
which is loaded into a `~jwst.datamodels.library.ModelLibrary` object,
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
     defined from a single exposure or from a dithered set of exposures
     can be found `at this link <https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-operations/nircam-dithers-and-mosaics>`_
   * Fill in pixels that have no valid contribution from any input exposure with the value
     specified by the ``fillval`` parameter.

#. If the ``save_intermediate_results`` parameter is set to `True`, write the resampled images to disk
   with the suffix ``_outlier_i2d.fits``. These are not saved if ``resample_data`` is set to `False` because
   they would be copies of the input models.

#. If resampling is turned off, use the input data itself in place of the resampled data
   for subsequent processing.

#. Construct and apply a bad pixel mask to the resampled data based on the drizzled weights.
   The ``weight_type`` parameter indicates the type of weighting image to apply with the bad pixel mask.
   The ``maskpt`` parameter sets the threshold weight value such that any pixel
   with a weight below this value gets flagged as bad.

#. Create a median image from all grouped observation mosaics pixel-by-pixel, i.e., averaging over groups.
   If ``save_intermediate_results`` is set to `True`, the median image is written out to disk with the
   suffix ``_median.fits``.

#. Blot (inverse of resampling) the median image back to match each original input image, and write
   the blotted images to disk with the suffix ``_blot.fits`` if ``save_intermediate_results`` is `True`.

#. Perform statistical comparison between blotted image and original image to identify outliers.

   * If resampling is disabled (``resample_data=False``), compare the median image directly
     to each input image and ignore all sub-bullets below this one.
     In this case, compute the outlier mask using the following formula:

     .. math::
        | image\_input - image\_median | > SNR \times input\_err

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

     .. math::
        | image\_input - image\_blotted | > scale \times image\_deriv + SNR \times noise

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
^^^^^^^^^^^^^^^^^^^^^
The outlier detection algorithm for imaging can require a lot of memory
depending on the number of inputs, the size of each input, and the size of the
final output product.  Specifically:

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
with the use of the ``in_memory`` parameter, which defaults to `True`.
The full impact of setting this parameter to `False` includes:

#. The input `~jwst.datamodels.library.ModelLibrary` object is loaded with ``on_disk=True``.
   This ensures that input models are loaded into memory one at at time,
   and saved to a temporary file when not in use; these read-write operations are handled internally by
   the `~jwst.datamodels.library.ModelLibrary` object.

#. Computing the median image works by writing the resampled data frames to appendable files
   on disk that are split into sections spatially but contain the entire ``groups``
   axis. The section size is set to use roughly the same amount of memory as a single resampled
   model, and since the resampled models are discarded from memory by the time the median calculation
   happens, this choice avoids increasing the overall memory usage of the step.
   Those sections are then read in one at a time to compute the median image.

These changes result in a minimum amount of memory usage during processing, but runtimes are
longer because many read and write operations are needed.
Note that if a `~jwst.datamodels.library.ModelLibrary` object
is input to the step, the memory behavior of the step is read from the ``on_disk`` status
of the `~jwst.datamodels.library.ModelLibrary` object,
and the ``in_memory`` parameter of the step is ignored.
When running :ref:`calwebb_image3 <calwebb_image3>`,
the ``in_memory`` flag should be set at the pipeline level
as the step-specific flag will be ignored, e.g.::

    strun calwebb_image3 asn.json --in-memory=False

.. _outlier-detection-spec:

Slit-like Spectroscopic Data
----------------------------

This section explains the interface for applying ``outlier_detection`` to slit-like
spectroscopic observations. The algorithm is very similar to the
:ref:`imaging algorithm <outlier-detection-imaging>`, and much of the same code is used.
A :ref:`Stage 3 association <asn-level3-techspecs>`,
which is loaded into a `~jwst.datamodels.container.ModelContainer` object,
serves as the input and output to this step, and the `~jwst.datamodels.container.ModelContainer`
is converted into a `~jwst.datamodels.library.ModelLibrary` object to allow sharing code
with the imaging mode.

This routine performs identical operations to the imaging mode, with the following *exceptions*:

#. Error thresholding is handled differently: The error arrays are resampled and median-combined
   along with the data arrays, and the median error image is used to identify outliers
   instead of the input error images for each exposure. This median error image is included
   alongside the median datamodel (in the ``err`` extension) if ``save_intermediate_results``
   is `True`.

#. Resampling is handled by a different class, `~jwst.resample.resample_spec.ResampleSpec`
   instead of `~jwst.resample.resample.ResampleImage`.

#. The resampled images are written out to disk with suffix "outlier_s2d" instead of
   "outlier_i2d" if the ``save_intermediate_results`` parameter is set to `True`.

#. The ``in_memory`` parameter has no effect, and all operations are performed in memory.

.. _outlier-detection-tso:

Time-Series Observations (TSO) Data
-----------------------------------

This section explains the interface for applying ``outlier_detection`` to time
series observations (TSO).

Normal imaging data benefit from combining all integrations into a
single image. TSO data's value, however, comes from looking for variations from one
integration to the next.  The outlier detection algorithm, therefore, gets run with
a few variations to accommodate the nature of these 3D data.
A `~stdatamodels.jwst.datamodels.CubeModel` object serves as the basic format for all
processing performed by this step. This routine performs the following operations:

#. Convert input data into a `~stdatamodels.jwst.datamodels.CubeModel` (3D data array)
   if a `~jwst.datamodels.container.ModelContainer` of 2D data arrays is provided.

#. Do not attempt resampling; data are assumed to be aligned and have an identical WCS.
   This is true automatically for a `~stdatamodels.jwst.datamodels.CubeModel`.

#. Apply a bad pixel mask to the input data based on the input DQ arrays and the ``good_bits``
   parameter.

#. Compute a median cube by combining all planes in the
   `~stdatamodels.jwst.datamodels.CubeModel` pixel-by-pixel using a
   rolling-median algorithm, in order to flag outliers integration-by-integration but
   preserve real time variability. The ``rolling_window_width`` parameter specifies the
   number of integrations over which to compute the median.

#. If the ``save_intermediate_results`` parameter is set to `True`, write the rolling-median
   `~stdatamodels.jwst.datamodels.CubeModel` to disk with the suffix ``_median.fits``.

#. Perform a statistical comparison frame-by-frame between the rolling-median cube and
   the input data. The formula used is the same as for imaging data without resampling:

   .. math::
      | image\_input - image\_median | > SNR \times input\_err

#. Update DQ arrays with flags and set SCI, ERR, and variance arrays to NaN at the location
   of identified outliers.

.. _outlier-detection-ifu:

Integral Field Unit (IFU) Data
------------------------------

This section explains the interface for applying ``outlier_detection`` to
Integral Field Unit (IFU) observations, like those taken with NIRSpec and MIRI.
A :ref:`Stage 3 association <asn-level3-techspecs>`,
which is loaded into a `~jwst.datamodels.container.ModelContainer` object,
serves as the basic format for all processing performed by this step.

After the JWST launch, it was discovered that the bad pixels on the MIRI detectors vary with time.
The pixels varied from usable to unusable, and at times, back to usable on a time scale that was too short
(sometimes as short as 2 days) to fold into the bad pixel mask applied in the
:ref:`calwebb_detector1 <calwebb_detector1>` pipeline. At this time it is believed that NIRSpec IFU data
also have bad pixels that vary with time, though the time variation is still under study.
The ``outlier_detection`` step is designed to flag these pixels as outliers, in addition
to cosmic ray hits that were not flagged by the :ref:`jump <jump_step>` step.

The basis of the outlier detection flagging for IFU data is to look for pixels on the detector
that are regularly discrepant from their neighbors, with a sharper division than could be explained
by the detector PSF.
This routine performs the following operations:

#. Extract parameter settings for the input `~jwst.datamodels.container.ModelContainer`
   and merge them with any user-provided values.
   See :ref:`outlier detection arguments <outlier_detection_step_args>` for the full list of parameters.

#. Loop over ``cal`` files, computing nearest-neighbor differences for each pixel
   in the along-dispersion direction.
   For MIRI, with the dispersion axis along the y axis, the neighbors that are used to
   to find the differences are to the left and right of each pixel being examined.
   For NIRSpec, with the dispersion along the x axis, the neighbors that are used to
   find the differences are above and below the pixel being examined.
   The smaller of the two (left/right or up/down) differences is stored as the difference value for each
   pixel; this avoids artifacts from bright edges.

#. Compare the nearest-neighbor differences across science exposures to find the minimum
   neighbor difference at each detector pixel.

#. Determine a local spatial median of the minimum difference array using a median filter with a kernel size
   set by the user according to the ``kernel_size`` parameter.

#. Normalize the minimum difference array by the local median.

#. Select outliers by flagging those normalized minimum values larger than the ``threshold_percent``
   parameter.

#. Update DQ arrays with flags and set SCI, ERR, and variance arrays to NaN at the location
   of identified outliers.

.. _outlier-detection-coron:

Coronagraphic Data
------------------

This section explains the interface for applying ``outlier_detection`` to coronagraphic
image observations. A `~stdatamodels.jwst.datamodels.CubeModel` serves as the basic format
for all processing performed by this step. This routine performs the following operations:

#. Extract parameter settings from input model and merge them with any user-provided values.
   See :ref:`outlier detection arguments <outlier_detection_step_args>` for the full list
   of parameters.

#. Do not attempt resampling; data are assumed to be aligned and have an identical WCS.
   This is true automatically for a `~stdatamodels.jwst.datamodels.CubeModel`.

#. Create a median image over the ``groups`` (exposures, planes of cube) axis,
   preserving the spatial (x, y) dimensions of the cube.

   * The ``maskpt`` parameter sets the percentage of the weight image values to
     use, and any pixel with a weight below this value gets flagged as "bad".

#. Perform statistical comparison between median image and original image to identify outliers.

   The core detection algorithm uses the following to generate an outlier mask:

   .. math::
      | image\_input - image\_median | > SNR \times input\_err

#. Update DQ arrays with flags and set SCI, ERR, and variance arrays to NaN at the location
   of identified outliers.

Reference Files
---------------
The ``outlier_detection`` step does not use any reference files.
