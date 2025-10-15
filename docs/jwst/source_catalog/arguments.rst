Arguments
=========

The ``source_catalog`` step uses the following user-settable arguments:

**General parameters:**

* ``--suffix``: A string value giving the file name suffix to use for
  the output catalog file. (Default='cat')

**Aperture photometry parameters:**

* ``--aperture_ee1``: An integer value of the smallest aperture
  encircled energy value. (Default=30)

* ``--aperture_ee2``: An integer value of the middle aperture encircled
  energy value. (Default=50)

* ``--aperture_ee3``: An integer value of the largest aperture encircled
  energy value. (Default=70)

* ``--ci1_star_threshold``: A floating-point value of the threshold
  compared to the concentration index calculated from aperture_ee1
  and aperture_ee2 that is used to determine whether a source is a
  star. Sources must meet the criteria of both ci1_star_threshold and
  ci2_star_threshold to be considered a star. (Default=2.0)

* ``--ci2_star_threshold``: A floating-point value of the threshold
  compared to the concentration index calculated from aperture_ee2
  and aperture_ee3 that is used to determine whether a source is a
  star. Sources must meet the criteria of both ci1_star_threshold and
  ci2_star_threshold to be considered a star. (Default=1.8)


**General source finding parameters:**

* ``starfinder``: A `str` indicating the source detection algorithm to use.
  Allowed values: ``'iraf'``, ``'dao'``, ``'segmentation'``. (Default= ``'segmentation'``)

* ``snr_threshold``: A `float` value indicating SNR threshold above the
  background. Required for all star finders. (Default=3.0)

* ``bkg_boxsize``: A positive `int` indicating the background mesh box size
  in pixels. (Default=1000)

* ``kernel_fwhm``: A `float` value indicating the Gaussian kernel FWHM in
  pixels. (Default=2.0)


**Additional source finding parameters for "segmentation":**

* ``npixels``: An `int` value indicating the minimum number of
  connected pixels that comprises a segment (Default=25)

* ``connectivity``: An `int` value indicating the connectivity defining the
  neighborhood of a pixel. Options are ``4``, i.e., connected pixels touch along edges,
  or ``8``, i.e, connected pixels touch along edges or corners (Default=8)

* ``nlevels``: An `int` value indicating the number of multi-thresholding
  levels for deblending (Default=32)

* ``contrast``: A `float` value indicating the fraction of total source flux
  an object must have to be deblended (Default=0.001)

* ``multithresh_mode``: A `str` indicating the multi-thresholding mode.
  Allowed values: ``'exponential'``, ``'linear'``, ``'sinh'``.
  (Default= ``'exponential'``)

* ``localbkg_width``: An `int` value indicating the width of rectangular
  annulus used to compute local background around each source. If set to 0,
  then local background will not be subtracted. (Default=0)

* ``apermask_method``: A `str` indicating the method used to handle
  neighboring sources when performing aperture photometry.
  Allowed values: ``'correct'``, ``'mask'``, ``'none'``. (Default= ``'correct'``)

* ``kron_params``: A tuple of `float` values indicating the
  parameters defining Kron aperture. If None,
  the parameters ``(2.5, 1.4, 0.0)`` are used. (Default=None)

* ``deblend``: A `bool` indicating whether to deblend sources. (Default=False)


**Additional source finding parameters for "dao" and "iraf":**

* ``minsep_fwhm``: A `float` value indicating the minimum separation between
  detected objects in units of number of FWHMs. (Default=0.0)

* ``sigma_radius``: A `float` value indicating the truncation radius of the
  Gaussian kernel in units of number of FWHMs. (Default=1.5)

* ``sharplo``: A `float` value indicating The lower bound on sharpness
  for object detection. (Default=0.5)

* ``sharphi``: A `float` value indicating the upper bound on sharpness
  for object detection. (Default=2.0)

* ``roundlo``: A `float` value indicating the lower bound on roundness
  for object detection. (Default=0.0)

* ``roundhi``: A `float` value indicating the upper bound on roundness
  for object detection. (Default=0.2)

* ``brightest``: A positive `int` value indicating the number of brightest
  objects to keep. If None, keep all objects above the threshold. (Default=200)

* ``peakmax``: A `float` value used to filter out objects with pixel values
  >= ``peakmax``. (Default=None)

.. warning::
  Different source finding algorithms have different values for the
  ``sharplo``, ``sharphi``, ``roundlo`` and ``roundhi`` parameters. These
  parameters should be adjusted to match the algorithm selected by the
  ``starfinder`` parameter. See documentation for
  `~photutils.detection.IRAFStarFinder`
  and `~photutils.detection.DAOStarFinder`.