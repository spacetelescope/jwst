Description
===========

Overview
--------
This step creates image catalogs of point-like sources whose
centroids are then used to compute corrections to the WCS of
the input images such that sky catalogs obtained from the
image catalogs using the corrected WCS will align on the sky.

Source Detection
----------------
Stars are detected in the image using the Photutils "daofind" function.
Photutils.daofind is an implementation of the `DAOFIND`_ algorithm
(`Stetson 1987, PASP 99, 191
<http://adsabs.harvard.edu/abs/1987PASP...99..191S>`_).  It searches
images for local density maxima that have a peak amplitude greater
than a specified threshold (the threshold is applied to a convolved
image) and have a size and shape similar to a defined 2D Gaussian
kernel.  Photutils.daofind also provides an estimate of the objects
roundness and sharpness, whose lower and upper bounds can be
specified.

.. _DAOFIND: http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?daofind

Alignment
---------
The source catalogs for each input image are compared to each other
and linear (affine) coordinate transformations that align these
catalogs are derived.

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
  be written out. (Default=`False`)

* ``catalog_format``: A `str` indicating catalog output file format.
  (Default='ecsv')

* ``kernel_fwhm``: A `float` value indicating the Gaussian kernel FWHM in
  pixels. (Default=2.5)

* ``snr_threshold``: A `float` value indicating SNR threshold above the
  background. (Default=5.0)

* ``brightest``: A positive `int` value indicating the number of brightest
  objects to keep. (Default=100)

* ``peakmax``: A `float` value used to filter out objects with pixel values
  >= ``peakmax``. (Default=None)

**Optimize alignment order:**

* ``enforce_user_order``: a boolean value indicating whether or not take the
  first image as a reference image and then align the rest of the images
  to that reference image in the order in which input images have been provided
  or to optimize order in which images are aligned. (Default=`False`)

**Reference Catalog parameters:**

* ``expand_refcat``: A boolean indicating whether or not to expand reference
  catalog with new sources from other input images that have been already
  aligned to the reference image. (Default=False)

**Object matching parameters:**

* ``minobj``: A positive `int` indicating minimum number of objects acceptable
  for matching. (Default=15)

* ``searchrad``: A `float` indicating the search radius in arcsec for a match.
  (Default=1.0)

* ``use2dhist``: A boolean indicating whether to use 2D histogram to find
  initial offset. (Default=True)

* ``separation``: Minimum object separation in arcsec. (Default=0.5)

* ``tolerance``: Matching tolerance for ``xyxymatch`` in arcsec. (Default=1.0)

* ``xoffset``: Initial guess for X offset in arcsec. (Default=0.0)

* ``yoffset``: Initial guess for Y offset in arcsec. (Default=0.0)

**Catalog fitting parameters:**

* ``fitgeometry``: A `str` value indicating the type of affine transformation
  to be considered when fitting catalogs. Allowed values: {``'shift'``,
  ``'rscale'``, ``'general'``}. (Default="general")

* ``nclip``: A non-negative integer number of clipping iterations
  to use in the fit. (Default = 3)

* ``sigma``: A positive `float` indicating the clipping limit, in sigma units,
  used when performing fit. (Default=3.0)

Further Documentation
---------------------
The underlying algorithms are described in more detail at

https://tweakwcs.readthedocs.io/en/latest/


Reference Files
===============
This step does not require any reference files.
