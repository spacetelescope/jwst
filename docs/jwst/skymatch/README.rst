Description
============


Overview
--------
The ``skymatch`` step can be used to compute sky values of input images or
it can be used to compute corrections that need to be applied to images such
as to "equalize" (match) sky background in input images.
When running ``skymatch`` step in a matching mode, ``skymatch`` compares
*total* signal levels in *the overlap regions*
(instead of doing this comparison on a per-pixel basis,
cf. :doc:`mrs_imatch step <../mrs_imatch/README>`) of a set of input images
and computes the signal offsets for each image that will minimize
the residuals across the entire set in the least squares sence. This comparison
is performed directly on input images without resampling them onto a common
grid. By default the sky value computed for each image is recorded, but
not actually subtracted from the images.


Assumptions
-----------

When matching sky background code needs to compute bounding polygon
intersections in world coordinates. Therefore, input images need to have
valid WCS.


Algorithm
---------
The ``skymatch`` step provides several methods for constant sky background
value computations.

First method, called ``'localmin'`` essentially is an enhanced version of the
original sky subtraction method used in older
`astrodrizzle <http://stsdas.stsci.edu/\
stsci_python_sphinxdocs_2.13/drizzlepac/astrodrizzle.html>`_\ versions. This
method simply computes the mean/median/mode/etc. value of the "sky" separately
in each input image. This method was upgraded to be able to use
DQ flags and user supplied masks to remove "bad" pixels from being
used for sky statistics computations. Values different from zero in
user-supplied masks indicate "good" data pixels.

In addition to the classical ``'localmin'``,
two other methods have been introduced: ``'globalmin'`` and
``'match'``, as well as a combination of the two -- ``'globalmin+match'``.

- The ``'globalmin'`` method computes the minimum sky value across *all*
  input images. That *single sky value* is then considered to be
  the background in *all input images*.

- The ``'match'`` algorithm computes constant (within an image) value
  corrections to be applied to input images such that the mismatch in computed
  background values between *all* pairs of images is minimized in the least
  squares sence. For each pair of images background mismatch is computed
  *only* in the regions in which the two images intersect.

  This makes ``'match'`` sky computation algorithm particularly useful
  for "equalizing" sky values in large mosaics in which one may have
  only (at least) pair-wise intersection of images without having
  a common intersection region (on the sky) in all images.

- The ``'globalmin+match'`` algorithm combines ``'match'`` and
  ``'globalmin'`` methods. It uses ``'globalmin'``
  algorithm to find a baseline sky value common to all input images
  and the ``'match'`` algorithm to "equalize" sky values among images.


Limitations and Discussions
---------------------------
Primary reason for introducing "sky match" algorithm was to try to
equalize the sky in large mosaics in which computation of the
"absolute" sky is difficult due to the presence of large diffuse
sources in the image. As discussed above, the skymatch step
accomplishes this by comparing "sky values" in input images in the
overlap regions (that is common to a pair of images). Quite obviously the
quality of sky "matching" will depend on how well these "sky values"
can be estimated. We use quotation marks around *sky values* because
for some image "true" background may not be present at all and the
measured sky may be the surface brightness of large galaxy, nebula, etc.

Here is a brief list of possible limitations/factors that can affect
the outcome of the matching (sky subtraction in general) algorithm:

* Since sky subtraction is performed on *flat-fielded* but
  *not distortion corrected* images, it is important to keep in mind
  that flat-fielding is performed to obtain uniform surface brightness
  and not flux. This distinction is important for images that have
  not been distortion corrected. As a consequence, it is advisable that
  point-like sources be masked through the user-supplied mask files.
  Values different from zero in user-supplied masks indicate "good" data
  pixels. Alternatively, one can use `upper` parameter to limit the use of
  bright objects in sky computations.

* Normally, distorted flat-fielded images contain cosmic rays. This
  algorithm does not perform CR cleaning. A possible way of minimizing
  the effect of the cosmic rays on sky computations is to use
  clipping (\ `nclip` > 0) and/or set `upper` parameter to a value
  larger than most of the sky background (or extended source) but
  lower than the values of most CR pixels.

* In general, clipping is a good way of eliminating "bad" pixels:
  pixels affected by CR, hot/dead pixels, etc. However, for
  images with complicated backgrounds (extended galaxies, nebulae,
  etc.), affected by CR and noise, clipping process may mask different
  pixels in different images. If variations in the background are
  too strong, clipping may converge to different sky values in
  different images even when factoring in the "true" difference
  in the sky background between the two images.

* In general images can have different "true" background values
  (we could measure it if images were not affected by large diffuse
  sources). However, arguments such as `lower` and `upper` will
  apply to all images regardless of the intrinsic differences
  in sky levels.


Reference Files
===============
This step does not require any reference files.

