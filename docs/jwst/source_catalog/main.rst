Description
===========

This step creates a final catalog of source photometry and morphologies.


Source Detection
^^^^^^^^^^^^^^^^

Sources are detected using `image segmentation
<https://en.wikipedia.org/wiki/Image_segmentation>`_, which is a
process of assigning a label to every pixel in an image such that
pixels with the same label are part of the same source.  The
segmentation procedure used is from `Photutils source extraction
<https://photutils.readthedocs.io/en/latest/segmentation.html>`_
and is called the threshold method, where detected sources must have a
minimum number of connected pixels that are each greater than a
specified threshold value in an image.  The threshold level is usually
defined at some multiple of the background standard deviation (sigma)
above the background.  The image can also be filtered before
thresholding to smooth the noise and maximize the detectability of
objects with a shape similar to the filter kernel.


Source Deblending
^^^^^^^^^^^^^^^^^

Note that overlapping sources are detected as single sources.
Separating those sources requires a deblending procedure, such as a
multi-thresholding technique used by `SExtractor
<https://www.astromatic.net/software/sextractor>`_.  Here we use the
`Photutils deblender
<https://photutils.readthedocs.io/en/latest/segmentation.html#source-deblending>`_,
which is an algorithm that deblends sources using a combination of
multi-thresholding and `watershed segmentation
<https://en.wikipedia.org/wiki/Watershed_(image_processing)>`_.  In
order to deblend sources, they must be separated enough such that
there is a saddle between them.


Source Photometry and Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After detecting sources using image segmentation, we can measure their
photometry, centroids, and morphological properties.  Here we use the
functions in `Photutils
<https://photutils.readthedocs.org/en/latest/segmentation.html>`_.  The
properties that are currently calculated for each source include
source centroids (both in pixel and sky coordinates), fluxes (and
errors), AB magnitudes (and errors), area, semimajor and semiminor
axis lengths, orientation, and sky coordinates at bounding box
corners.
