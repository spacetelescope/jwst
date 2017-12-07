
Description
===========

This step creates a catalog of point-like sources whose centroids are
then used to aligning images using the Tweakreg algorithm.

Source detection
----------------

Stars are detected in the image using Photutils' "daofind" function.
Photutils.daofind is an implementation of the `DAOFIND`_ algorithm
(`Stetson 1987, PASP 99, 191
<http://adsabs.harvard.edu/abs/1987PASP...99..191S>`_).  It searches
images for local density maxima that have a peak amplitude greater
than a specified threshold (the threshold is applied to a convolved
image) and have a size and shape similar to a defined 2D Gaussian
kernel.  Photutils.daofind also provides an estimate of the objects'
roundness and sharpness, whose lower and upper bounds can be
specified.

.. _DAOFIND: http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?daofind

Alignment
---------

The source catalog for each input image is used to compute the affine
transformation that aligns the images to each other.  This affine
transformation is appended to the GWCS object of each image.