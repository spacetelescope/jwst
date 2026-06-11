Step Arguments
==============

The ``source_catalog`` step has the following optional arguments.

**General parameters:**

``--suffix`` (str, default='cat')
  File name suffix to use for the output catalog file.

**Aperture photometry parameters:**

``--aperture_ee1`` (int, default=30)
  Value of the smallest aperture encircled energy.

``--aperture_ee2`` (int, default=50)
  Value of the middle aperture encircled energy.

``--aperture_ee3`` (int, default=70)
  Value of the largest aperture encircled energy.

``--ci1_star_threshold`` (float, default=2.0)
  The threshold
  compared to the concentration index calculated from ``aperture_ee1``
  and ``aperture_ee2`` that is used to determine whether a source is a
  star. Sources must meet the criteria of both ``ci1_star_threshold`` and
  ``ci2_star_threshold`` to be considered a star.

``--ci2_star_threshold`` (float, default=1.8)
  The threshold
  compared to the concentration index calculated from ``aperture_ee2``
  and ``aperture_ee3`` that is used to determine whether a source is a
  star. Sources must meet the criteria of both ``ci1_star_threshold`` and
  ``ci2_star_threshold`` to be considered a star.

**General source finding parameters:**

``--starfinder`` (str, default='segmentation')
  The source detection algorithm to use.
  Allowed values: ``'iraf'``, ``'dao'``, ``'segmentation'``.

``--snr_threshold`` (float, default=3.0)
  SNR threshold above the
  background. Required for all star finders.

``--bkg_boxsize`` (int, default=1000)
  A positive value indicating the background mesh box size
  in pixels.

``--kernel_fwhm`` (float, default=2.0)
  The Gaussian kernel FWHM in pixels.

**Additional source finding parameters for "segmentation":**

``--npixels`` (int, default=25)
  The minimum number of connected pixels that comprises a segment.

``--connectivity`` (int, default=8)
  The connectivity defining the neighborhood of a pixel.
  Options are ``4``, i.e., connected pixels touch along edges,
  or ``8``, i.e, connected pixels touch along edges or corners.

``--nlevels`` (int, default=32)
  The number of multi-thresholding levels for deblending.

``--contrast`` (float, default=0.001)
  The fraction of total source flux an object must have to be deblended.

``--multithresh_mode`` (str, default='exponential')
  The multi-thresholding mode.
  Allowed values: ``'exponential'``, ``'linear'``, ``'sinh'``.

``--localbkg_width`` (int, default=0)
  The width of rectangular
  annulus used to compute local background around each source. If set to 0,
  then local background will not be subtracted.

``--apermask_method`` (str, default='correct')
  The method used to handle
  neighboring sources when performing aperture photometry.
  Allowed values: ``'correct'``, ``'mask'``, ``'none'``.

``--kron_params`` (tuple of float, default=None)
  The parameters defining Kron aperture. If None,
  the parameters ``(2.5, 1.4, 0.0)`` are used.

``--deblend`` (bool, default=True)
  Indicator for whether to deblend sources.

**Additional source finding parameters for "dao" and "iraf":**

``--minsep_fwhm`` (float, default=0.0)
  The minimum separation between
  detected objects in units of number of FWHMs.

``--sigma_radius`` (float, default=1.5)
  The truncation radius of the
  Gaussian kernel in units of number of FWHMs.

``--sharplo`` (float, default=0.5)
  The lower bound on sharpness for object detection.

``--sharphi`` (float, default=2.0)
  The upper bound on sharpness for object detection.

``--roundlo`` (float, default=0.0)
  The lower bound on roundness for object detection.

``--roundhi`` (float, default=0.2)
  The upper bound on roundness for object detection.

``--brightest`` (int, default=200)
  A positive value indicating the number of brightest
  objects to keep. If None, keep all objects above the threshold.

``--peakmax`` (float, default=None)
  A value used to filter out objects with pixel values
  ``>= peakmax``.

.. warning::
  Different source finding algorithms have different values for the
  ``sharplo``, ``sharphi``, ``roundlo``, and ``roundhi`` parameters. These
  parameters should be adjusted to match the algorithm selected by the
  ``starfinder`` parameter. See documentation for
  `~photutils.detection.IRAFStarFinder`
  and `~photutils.detection.DAOStarFinder`.
