Description
===========

:Class: `jwst.tweakreg.TweakRegStep`
:Alias: tweakreg

Overview
--------
This step creates image catalogs of point-like sources whose
centroids are then used to compute corrections to the WCS of
the input images such that sky catalogs obtained from the
image catalogs using the corrected WCS will align on the sky.

Source Detection
----------------
If the ``meta.tweakreg_catalog`` attribute of input data models is a non-empty
string and ``use_custom_catalogs`` is `True`, then it will be interpreted
as a file name of a user-provided source catalog. The catalog must be in a
format automatically recognized by :py:meth:`~astropy.table.Table.read`.

When the ``meta.tweakreg_catalog`` attribute of input data models is `None` or
an empty string, then the ``tweakreg`` step will attempt to detect sources in the
input images. Stars are detected in the image with one of the following source
detection algorithms: ``photutils.detection.DAOStarFinder`` (default),
``photutils.detection.IRAFStarFinder``, or ``photutils.segmentation.SourceFinder``
in conjunction with ``photutils.segmentation.SourceCatalog``.

DAOStarFinder is an implementation of the `DAOFIND`_ algorithm
(`Stetson 1987, PASP 99, 191
<http://adsabs.harvard.edu/abs/1987PASP...99..191S>`_).  It searches
images for local density maxima that have a peak amplitude greater
than a specified threshold (the threshold is applied to a convolved
image) and have a size and shape similar to a defined 2D Gaussian
kernel.  DAOFind also provides an estimate of the object's
roundness and sharpness, whose lower and upper bounds can be
specified.

IRAFStarFinder is a Python implementation of the IRAF star finding algorithm,
which also calculates the objects' centroids, roundness, and sharpness.
However, IRAFStarFinder uses image moments
instead of 1-D Gaussian fits to projected light distributions like
DAOStarFinder.

SourceFinder implements a segmentation algorithm that identifies
sources in an image based on a number of connected pixels above a
specified threshold value.  The sources are deblended using a
combination of multi-thresholding and watershed segmentation.
SourceCatalog finds the centroids of these sources, which are used
as the retrieved star positions.

.. warning::
    It has been shown (`STScI Technical Report JWST-STScI-008116, SM-12
    <https://www.stsci.edu/~goudfroo/NIRISSdoc/Centroid_Accuracies_Precisions_NIRISS_v2.pdf>`_)
    that for undersampled PSFs, e.g. for short-wavelength NIRISS
    imaging data, ``DAOStarFinder`` gives bad results no matter the input parameters
    due to its use of 1-D Gaussian fits.
    ``IRAFStarFinder`` or ``SourceFinder`` should be used instead.

.. note::
    ``SourceFinder`` is likely to detect non-stellar sources
    such as galaxies because sources are not assumed to be
    point-source-like. This may lead to mismatches between the
    derived source catalog and the reference catalog during the
    alignment step.

.. _DAOFIND: http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?daofind

Custom Source Catalogs
----------------------
Source detection built into the ``tweakreg`` step can be disabled by
providing a file name to a custom source catalog in the
``meta.tweakreg_catalog`` attribute of input data models.
The catalog must be in a format automatically recognized by
:py:meth:`~astropy.table.Table.read`. The catalog must contain
either ``'x'`` and ``'y'`` or ``'xcentroid'`` and ``'ycentroid'`` columns which
indicate source *image* coordinates (in pixels). Pixel coordinates are
0-indexed. An optional column in the catalog is the ``'weight'`` column,
which when present, will be used in fitting.

For the ``tweakreg`` step to use user-provided input source catalogs,
``use_custom_catalogs`` parameter of the ``tweakreg`` step must be set to
`True`.

In addition to setting the ``meta.tweakreg_catalog`` attribute of input data
models to the custom catalog file name, the ``tweakreg_step`` also supports two
other ways of supplying custom source catalogs to the step:

1. Adding ``tweakreg_catalog`` attribute to the ``members`` of the input ASN
   table - see `~jwst.datamodels.ModelLibrary` for more details.
   Catalog file names are relative to ASN file path.

2. Providing a simple two-column text file, specified via step's parameter
   ``catfile``, that contains input data models' file names in the first column
   and the file names of the corresponding catalogs in the second column.
   Catalog file names are relative to ``catfile`` file path.

Specifying custom source catalogs via either the input ASN file or
``catfile`` will update input data models' ``meta.tweakreg_catalog``
attributes to the catalog file names provided in either in the ASN file or
``catfile``.

.. note::
    When custom source catalogs are provided via both ``catfile`` and
    ASN file members' attributes, the ``catfile`` takes precedence and
    catalogs specified via ASN file are ignored altogether.

.. note::
    1. Providing a data model file name in the ``catfile`` and leaving
       the corresponding source catalog file name empty -- same as setting
       ``'tweakreg_catalog'`` in the ASN file to an empty string ``""`` --
       would set the corresponding input data model's ``meta.tweakreg_catalog``
       attribute to `None`. In this case, ``tweakreg_step`` will automatically
       generate a source catalog for that data model.

    2. If an input data model is not listed in the ``catfile`` or does not
       have the ``'tweakreg_catalog'`` attribute provided in the ASN file,
       then the catalog file name in that model's ``meta.tweakreg_catalog``
       attribute will be used. If ``model.meta.tweakreg_catalog`` is `None`,
       ``tweakreg_step`` will automatically generate a source catalog for
       that data model.

Alignment
---------
The source catalogs for each input image are compared to each other
and linear (affine) coordinate transformations that align these
catalogs are derived.  This fit ensures that all the input images
are aligned relative to each other.  This step produces a combined
source catalog for the entire set of input images as if they were
combined into a single mosaic.

If the step parameter ``abs_refcat`` is set to 'GAIADR3', 'GAIADR2', or 'GAIADR1',
an astrometric reference catalog then gets generated by querying
a GAIA-based astrometric catalog web service for all astrometrically
measured sources in the combined field-of-view of the set of input
images. This catalog is generated from the catalogs available
through the `STScI MAST Catalogs`_ and has the ability to account
for proper motion to a given epoch. The epoch is computed from the observation date and time
of the input data. If ``abs_refcat`` is set to a path to an existing
file, i.e., a user-supplied external reference catalog,
then the catalog will be read from that file. The catalog must be readable
into an :py:meth:`~astropy.table.Table` object and contain either
``'RA'`` and ``'DEC'`` columns (in degrees) or an Astropy-readable ``sky_centroid``.
An optional column in the catalog is the ``'weight'`` column, which when present,
will be used in fitting.

.. _STScI MAST Catalogs: https://outerspace.stsci.edu/display/MASTDATA/Catalog+Access

The combined source catalog derived in the first step
then gets cross-matched and fit to this astrometric reference catalog.
The results of this one fit then gets back-propagated to all the
input images to align them all to the astrometric reference frame while
maintaining the relative alignment between the images.


Grouping
--------

Images taken at the same time (e.g., NIRCam images from all short-wave
detectors) can be aligned together; that is, a single correction
can be computed and applied to all these images because any error in
telescope pointing will be identical in all these images and it is assumed
that the relative positions of (e.g., NIRCam) detectors do not change.
Identification of images that belong to the same "exposure" and therefore
can be grouped together is based on several attributes described in
`~jwst.datamodels.ModelLibrary`. This grouping is performed automatically
in the ``tweakreg`` step using the
`~jwst.datamodels.ModelLibrary.group_names` property.


However, when detector calibrations are not accurate, alignment of groups
of images may fail (or result in poor alignment). In this case, it may be
desirable to align each image independently. This can be achieved either by
setting the ``image_model.meta.group_id`` attribute to a unique string or integer
value for each image, or by adding the ``group_id`` attribute to the ``members`` of the input ASN
table - see `~jwst.datamodels.ModelLibrary` for more details.

.. note::
    Group ID (``group_id``) is used by both ``tweakreg`` and ``skymatch`` steps
    and so modifying it for one step will affect the results in another step.
    If it is desirable to apply different grouping strategies to the ``tweakreg``
    and ``skymatch`` steps, one may need to run each step individually and
    provide a different ASN as input to each step.

WCS Correction
--------------
The linear coordinate transformation computed in the previous step
is used to define tangent-plane corrections that need to be applied
to the GWCS pipeline in order to correct input image WCS.
This correction is implemented by inserting a ``v2v3corr`` frame with
tangent plane corrections into the GWCS pipeline of the image's WCS.

Step Arguments
--------------
The ``tweakreg`` step has the following optional arguments:

**Source finding parameters:**

* ``save_catalogs``: A boolean indicating whether or not the catalogs should
  be written out. This parameter is ignored for input data models whose
  ``meta.tweakreg_catalog`` is a non-empty string pointing to a user-supplied
  source catalog. (Default=False)

* ``use_custom_catalogs``: A boolean that indicates whether
  to ignore source catalog in the input data model's ``meta.tweakreg_catalog``
  attribute. If `False`, new catalogs will be generated by the ``tweakreg``
  step. (Default=False)

* ``catalog_format``: A `str` indicating catalog output file format.
  (Default= `'ecsv'`)

* ``catfile``: Name of the file with a list of custom user-provided catalogs.
  The file must contain a two-column list of format
  ``<input file name> <catalog file name>`` with one entry per input filename
  in the input association.
  This parameter has no effect if ``use_custom_catalogs`` is `False`.
  (Default= `''`)

* ``bkg_boxsize``: A positive `int` indicating the background mesh box size
  in pixels. (Default=400)

* ``starfinder``: A `str` indicating the source detection algorithm to use.
  Allowed values: `'iraf'`, `'dao'`, `'segmentation'`. (Default= `'iraf'`)

* ``snr_threshold``: A `float` value indicating SNR threshold above the
  background. Required for all star finders. (Default=10.0)

* ``kernel_fwhm``: A `float` value indicating the Gaussian kernel FWHM in
  pixels. (Default=2.5)

**Additional source finding parameters for DAO and IRAF:**

* ``minsep_fwhm``: A `float` value indicating the minimum separation between
  detected objects in units of number of FWHMs. (Default=0.0)

* ``sigma_radius``: A `float` value indicating the truncation radius of the
  Gaussian kernel in units of number of FWHMs. (Default=2.5)

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
  [IRAFStarFinder](https://photutils.readthedocs.io/en/stable/api/photutils.detection.IRAFStarFinder.html)
  and [DAOStarFinder](https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html).

**Additional source finding parameters for segmentation:**

* ``npixels``: An `int` value indicating the minimum number of
  connected pixels that comprises a segment (Default=10)

* ``connectivity``: An `int` value indicating the connectivity defining the
  neighborhood of a pixel. Options are `4`, i.e., connected pixels touch along edges,
  or `8`, i.e, connected pixels touch along edges or corners (Default=8)

* ``nlevels``: An `int` value indicating the number of multi-thresholding
  levels for deblending (Default=32)

* ``contrast``: A `float` value indicating the fraction of total source flux
  an object must have to be deblended (Default=0.001)

* ``multithresh_mode``: A `str` indicating the multi-thresholding mode.
  Allowed values: `'exponential'`, `'linear'`, `'sinh'`.
  (Default= `'exponential'`)

* ``localbkg_width``: An `int` value indicating the width of rectangular
  annulus used to compute local background around each source. If set to 0,
  then local background will not be subtracted. (Default=0)

* ``apermask_method``: A `str` indicating the method used to handle
  neighboring sources when performing aperture photometry.
  Allowed values: `'correct'`, `'mask'`, `'none'`. (Default= `'correct'`)

* ``kron_params``: A tuple of `float` values indicating the
  parameters defining Kron aperture. If None,
  the parameters `(2.5, 1.4, 0.0)` are used. (Default=None)

**Optimize alignment order:**

* ``enforce_user_order``: a boolean value indicating whether or not take the
  first image as a reference image and then align the rest of the images
  to that reference image in the order in which input images have been provided
  or to optimize order in which images are aligned. (Default=False)

**Reference Catalog parameters:**

* ``expand_refcat``: A boolean indicating whether or not to expand reference
  catalog with new sources from other input images that have been already
  aligned to the reference image. (Default=False)

**Object matching parameters:**

* ``minobj``: A positive `int` indicating minimum number of objects acceptable
  for matching. (Default=15)

* ``searchrad``: A `float` indicating the search radius in arcsec for a match.
  (Default=2.0)

* ``use2dhist``: A boolean indicating whether to use 2D histogram to find
  initial offset. (Default=True)

* ``separation``: Minimum object separation in arcsec. It **must be** at least
  ``sqrt(2)`` times larger than ``tolerance``. (Default=1.0)

* ``tolerance``: Matching tolerance for ``xyxymatch`` in arcsec. (Default=0.7)

* ``xoffset``: Initial guess for X offset in arcsec. (Default=0.0)

* ``yoffset``: Initial guess for Y offset in arcsec. (Default=0.0)

**Catalog fitting parameters:**

* ``fitgeometry``: A `str` value indicating the type of affine transformation
  to be considered when fitting catalogs. Allowed values:

  - ``'shift'``: x/y shifts only
  - ``'rshift'``: rotation and shifts
  - ``'rscale'``: rotation and scale
  - ``'general'``: shift, rotation, and scale

  The default value is "rshift".

  .. note::
      Mathematically, alignment of images observed in different tangent planes
      requires ``fitgeometry='general'`` in order to fit source catalogs
      in the different images even if misalignment is caused only by a shift
      or rotation in the tangent plane of one of the images.

      However, under certain circumstances, such as small alignment errors or
      minimal dithering during observations that keep tangent planes of the
      images to be aligned almost parallel, then it may be more robust to
      use a ``fitgeometry`` setting with fewer degrees of freedom such as
      ``'rshift'``, especially for "ill-conditioned" source catalogs such as
      catalogs with very few sources, or large errors in source positions, or
      sources placed along a line or bunched in a corner of the image (not
      spread across/covering the entire image).

* ``nclip``: A non-negative integer number of clipping iterations
  to use in the fit. (Default=3)

* ``sigma``: A positive `float` indicating the clipping limit, in sigma units,
  used when performing fit. (Default=3.0)

**Absolute Astrometric fitting parameters:**

Parameters used for absolute astrometry to a reference catalog.

* ``abs_refcat``: String indicating what astrometric catalog should be used.
  Currently supported options: 'GAIADR1', 'GAIADR2', 'GAIADR3', a path to an existing
  reference catalog, `None`, or `''`. See
  :py:data:`jwst.tweakreg.tweakreg_step.SINGLE_GROUP_REFCAT`
  for an up-to-date list of supported built-in reference catalogs.

  When ``abs_refcat`` is `None` or an empty string, alignment to the
  absolute astrometry catalog will be turned off.
  (Default= `''`)

* ``abs_minobj``: A positive `int` indicating minimum number of objects
  acceptable for matching. (Default=15)

* ``abs_searchrad``: A `float` indicating the search radius in arcsec for
  a match. It is recommended that a value larger than ``searchrad`` be used for
  this parameter (e.g. 3 times larger) (Default=6.0)

* ``abs_use2dhist``: A boolean indicating whether to use 2D histogram to find
  initial offset. It is strongly recommended setting this parameter to `True`.
  Otherwise the initial guess for the offsets will be set to zero
  (Default=True)

* ``abs_separation``: Minimum object separation in arcsec. It **must be** at
  least ``sqrt(2)`` times larger than ``abs_tolerance``. (Default=1.0)

* ``abs_tolerance``: Matching tolerance for ``xyxymatch`` in arcsec.
  (Default=0.7)

* ``abs_fitgeometry``: A `str` value indicating the type of affine
  transformation to be considered when fitting catalogs. Allowed values:

  - ``'shift'``: x/y shifts only
  - ``'rshift'``: rotation and shifts
  - ``'rscale'``: rotation and scale
  - ``'general'``: shift, rotation, and scale

  The default value is "rshift". Note that the same conditions/restrictions
  that apply to ``fitgeometry`` also apply to ``abs_fitgeometry``.

* ``abs_nclip``: A non-negative integer number of clipping iterations
  to use in the fit. (Default=3)

* ``abs_sigma``: A positive `float` indicating the clipping limit, in sigma
  units, used when performing fit. (Default=3.0)

* ``save_abs_catalog``: A boolean specifying whether or not to write out the
  astrometric catalog used for the fit as a separate product. (Default=False)

**SIP approximation parameters:**

Parameters used to provide a SIP-based approximation to the WCS,
for FITS display. These parameter values should match the ones used
in the ``assign_wcs`` step.

* ``sip_approx``: A boolean flag to enable the computation of a SIP
  approximation. (Default=True)

* ``sip_degree``: A positive `int`, specifying the polynomial degree for
  the forward SIP fit. `None` uses the best fit; the maximum value allowed
  is 6. (Default=None)

* ``sip_max_pix_error``: A positive `float`, specifying the maximum
  error for the SIP forward fit, in units of pixels. Ignored if
  ``sip_degree`` is set to an explicit value. (Default=0.01)

* ``sip_inv_degree``: A positive `int`, specifying the polynomial degree for
  the inverse SIP fit. `None` uses the best fit; the maximum value allowed
  is 6. (Default=None)

* ``sip_max_inv_pix_error``: A positive `float`, specifying the maximum
  error for the SIP inverse fit, in units of pixels. Ignored if
  ``sip_inv_degree`` is set to an explicit value. (Default=0.01)

* ``sip_npoints``: Number of points for the SIP fit. (Default=12).

**stpipe general options:**

* ``output_use_model``: A boolean indicating whether to use `DataModel.meta.filename`
  when saving the results. (Default=True)

* ``in_memory``: A boolean indicating whether to keep models in memory, or to save
  temporary files on disk while not in use to save memory. (Default=True)

Further Documentation
---------------------
The underlying algorithms as well as formats of source catalogs are described
in more detail at

https://tweakwcs.readthedocs.io/en/latest/

Further description of the input parameters and algorithms for star finding
can be found at the following links:

* `DAOStarFinder`_
* `IRAFStarFinder`_
* `SourceFinder`_
* `SourceCatalog`_

.. _DAOStarFinder: https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html
.. _IRAFStarFinder: https://photutils.readthedocs.io/en/stable/api/photutils.detection.IRAFStarFinder.html
.. _SourceFinder: https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.SourceFinder.html
.. _SourceCatalog: https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.SourceCatalog.html

The alignment and WCS correction portions of the step are handled by the `stcal` package.
Additional documentation may be found `here <https://stcal.readthedocs.io/en/latest/stcal/tweakreg/index.html>`_.


Reference Files
===============
The ``tweakreg`` step uses the PARS-TWEAKREGSTEP parameter reference file.

.. include:: ../references_general/pars-tweakregstep_reffile.inc
