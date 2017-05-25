
Description
===========

This step creates a catalog of point-like sources whose centroids are
used to aid in aligning images in the "tweakreg" step.

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
