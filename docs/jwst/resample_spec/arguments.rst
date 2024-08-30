.. _resample_spec_step_args:

Step Arguments
==============
The `resample_spec` step has the following optional arguments that control
the behavior of the processing and the characteristics of the resampled
image.

``--pixfrac`` (float, default=1.0)
    The fraction by which input pixels are "shrunk" before being drizzled
    onto the output image grid, given as a real number between 0 and 1.

``--kernel`` (str, default='square')
    The form of the kernel function used to distribute flux onto the output
    image.  Available kernels for spectral data are `square` and `point`.
    Other drizzle kernels, available for imaging data, do not conserve
    spectral flux.

``--pixel_scale_ratio`` (float, default=1.0)
    Ratio of input to output spatial pixel scale.

    Values greater than 1 indicate that the input pixels have a larger spatial
    scale, so more output pixels will sample the same input pixel.  For example,
    a value of 2.0 means the output image would have 2 pixels sampling each input
    spatial pixel. If the input data has units of flux density (MJy/pixel),
    the output flux per pixel will be half the input flux per pixel.
    If the input data has units of surface brightness (MJy/sr), the output
    flux per pixel is not scaled.

    Note that this parameter is only applied in the cross-dispersion
    direction: sampling wavelengths are not affected.

    Ignored when ``pixel_scale`` or ``output_wcs`` are provided.

    .. note::
        If this parameter is modified from the default value, the extraction
        aperture for the :ref:`extract_1d <extract_1d_step>` step must
        also be modified, since it is specified in pixels.

``--pixel_scale`` (float, default=None)
    Absolute pixel scale in ``arcsec``. When provided, overrides
    ``pixel_scale_ratio``. Ignored when ``output_wcs`` is provided.

    If the input data has units of flux density (MJy/pixel), the output flux per
    pixel will be scaled by the ratio of the selected output pixel scale to an average
    input pixel scale.
    If the input data has units of surface brightness (MJy/sr),
    the output flux per pixel is not scaled.

    Note that this parameter is only applied in the cross-dispersion
    direction: sampling wavelengths are not affected.

    .. note::
        If this parameter is modified from the default value, the extraction
        aperture for the :ref:`extract_1d <extract_1d_step>` step must
        also be modified, since it is specified in pixels.

``--output_shape`` (tuple of int, default=None)
    Shape of the image (data array) using "standard" ``nx`` first and ``ny``
    second (opposite to the ``numpy.ndarray`` convention - ``ny`` first and
    ``nx`` second). This value will be assigned to
    ``pixel_shape`` and ``array_shape`` properties of the returned
    WCS object. When supplied from command line, it should be a comma-separated
    list of integers ``nx, ny``.

    .. note::
        Specifying ``output_shape`` *is required* when the WCS in
        ``output_wcs`` does not have ``bounding_box`` property set.

``--output_wcs`` (str, default='')
    File name of an ``ASDF`` file with a GWCS stored under the ``"wcs"`` key
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
        When ``output_wcs`` is specified, WCS-related arguments
        ``pixel_scale_ratio`` and ``pixel_scale`` will be ignored.

``--fillval`` (str, default='NAN')
    The value to assign to output pixels that have zero weight or do not
    receive any flux from any input pixels during drizzling.

``--weight_type`` (str, default='ivm')
    The weighting type for each input image.
    If `weight_type=ivm` (the default), the scaling value
    will be determined per-pixel using the inverse of the read noise
    (VAR_RNOISE) array stored in each input image. If the VAR_RNOISE array does
    not exist, the variance is set to 1 for all pixels (equal weighting).
    If `weight_type=exptime`, the scaling value will be set equal to the
    measurement time (TMEASURE) found in the image header if available;
    if unavailable, the scaling will be set equal to the exposure time (EFFEXPTM).

``--single`` (bool, default=False)
    If set to `True`, resample each input image into a separate output.  If
    `False` (the default), each input is resampled additively (with weights) to
    a common output.

``--blendheaders`` (bool, default=True)
    Blend metadata from all input images into the resampled output image.

``--in_memory`` (boolean, default=True)
  Specifies whether or not to load and create all images that are used during
  processing into memory. If ``False``, input files are loaded from disk when
  needed and all intermediate files are stored on disk, rather than in memory.
