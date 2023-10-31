.. _resample_step_args:

Step Arguments
==============
The `resample` step has the following optional arguments that control
the behavior of the processing and the characteristics of the resampled
image.

``--pixfrac`` (float, default=1.0)
    The fraction by which input pixels are "shrunk" before being drizzled
    onto the output image grid, given as a real number between 0 and 1.

``--kernel`` (str, default='square')
    The form of the kernel function used to distribute flux onto the output
    image.  Available kernels are `square`, `gaussian`, `point`, `tophat`,
    `turbo`, `lanczos2`, and `lanczos3`.

``--pixel_scale_ratio`` (float, default=1.0)
    Ratio of input to output pixel scale.  A value of 0.5 means the output
    image would have 4 pixels sampling each input pixel.
    Ignored when ``pixel_scale`` or ``output_wcs`` are provided.

``--pixel_scale`` (float, default=None)
    Absolute pixel scale in ``arcsec``. When provided, overrides
    ``pixel_scale_ratio``. Ignored when ``output_wcs`` is provided.

``--rotation`` (float, default=None)
    Position angle of output image’s Y-axis relative to North.
    A value of 0.0 would orient the final output image to be North up.
    The default of `None` specifies that the images will not be rotated,
    but will instead be resampled in the default orientation for the camera
    with the x and y axes of the resampled image corresponding
    approximately to the detector axes. Ignored when ``pixel_scale``
    or ``output_wcs`` are provided.

``--crpix`` (tuple of float, default=None)
    Position of the reference pixel in the image array in the ``x, y`` order.
    If ``crpix`` is not specified, it will be set to the center of the bounding
    box of the returned WCS object. When supplied from command line, it should
    be a comma-separated list of floats. Ignored when ``output_wcs``
    is provided.

``--crval`` (tuple of float, default=None)
    Right ascension and declination of the reference pixel. Automatically
    computed if not provided. When supplied from command line, it should be a
    comma-separated list of floats. Ignored when ``output_wcs`` is provided.

``--output_shape`` (tuple of int, default=None)
    Shape of the image (data array) using "standard" ``nx`` first and ``ny``
    second (as opposite to the ``numpy.ndarray`` convention - ``ny`` first and
    ``nx`` second). This value will be assigned to
    ``pixel_shape`` and ``array_shape`` properties of the returned
    WCS object. When supplied from command line, it should be a comma-separated
    list of integers ``nx, ny``.

    .. note::
        Specifying ``output_shape`` *is required* when the WCS in
        ``output_wcs`` does not have ``bounding_box`` property set.

``--output_wcs`` (str, default='')
    File name of a ``ASDF`` file with a GWCS stored under the ``"wcs"`` key
    under the root of the file. The output image size is determined from the
    bounding box of the WCS (if any). Argument ``output_shape`` overrides
    computed image size and it is required when output WCS does not have
    ``bounding_box`` property set or if ``pixel_shape`` or ``array_shape`` keys
    (see below) are not provided.

    Additional information may be stored under
    other keys under the root of the file. Currently, the following keys are
    recognized:

    - ``pixel_area``: Indicates average pixel area of the output WCS in
      units of steradians. When provided, this value will be used for updating
      photometric quantities  ``PIXAR_SR`` and ``PIXAR_A2`` of the output image.
      If ``pixel_area`` is not provided, the code will attempt to estimate
      this value from the WCS.

    - ``pixel_shape``: dimensions of the output image in the order (nx, ny).
      Overrides the value of ``array_shape`` if provided.

    - ``array_shape``: shape of the output image in ``numpy`` order: (ny, nx).

    .. note::
        When ``output_wcs`` is specified, WCS-related arguments such as
        ``pixel_scale_ratio``, ``pixel_scale``, ``rotation``, ``crpix``,
        and ``crval`` will be ignored.

``--fillval`` (str, default='INDEF')
    The value to assign to output pixels that have zero weight or do not
    receive any flux from any input pixels during drizzling.

``--weight_type`` (str, default='ivm')
    The weighting type for each input image.
    If `weight_type=ivm` (the default), the scaling value
    will be determined per-pixel using the inverse of the read noise
    (VAR_RNOISE) array stored in each input image. If the VAR_RNOISE array does
    not exist, the variance is set to 1 for all pixels (equal weighting).
    If `weight_type=exptime`, the scaling value will be set equal to the
    exposure time found in the image header.

``--single`` (bool, default=False)
    If set to `True`, resample each input image into a separate output.  If
    `False` (the default), each input is resampled additively (with weights) to
    a common output

``--blendheaders`` (bool, default=True)
    Blend metadata from all input images into the resampled output image.

``--allowed_memory`` (float, default=None)
    Specifies the fractional amount of free memory to allow when creating the
    resampled image. If ``None``, the environment variable
    ``DMODEL_ALLOWED_MEMORY`` is used. If not defined, no check is made. If the
    resampled image would be larger than specified, an ``OutputTooLargeError``
    exception will be generated.

    For example, if set to ``0.5``, only resampled images that use less than
    half the available memory can be created.
