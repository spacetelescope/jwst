Description
===========

This step creates a final catalog of source photometry and morphologies.


Source Detection
^^^^^^^^^^^^^^^^

Sources are detected using `image segmentation
<http://en.wikipedia.org/wiki/Image_segmentation>`_, which is a
process of assigning a label to every pixel in an image such that
pixels with the same label are part of the same source.  The
segmentation procedure used is from `Photutils
<http://photutils.readthedocs.org/en/latest/photutils/detection.html#source-extraction-using-image-segmentation>`_
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
<http://www.astromatic.net/software/sextractor>`_.  Here we use the
Photutils deblender, which is an experimental algorithm that deblends
sources using a combination of multi-thresholding and `watershed
segmentation
<https://en.wikipedia.org/wiki/Watershed_(image_processing)>`_.  In
order to deblend sources, they must be separated enough such that
there is a saddle between them.


Source Photometry and Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After detecting sources using image segmentation, we can measure their
photometry, centroids, and morphological properties.  Here we use the
functions in `Photutils
<http://photutils.readthedocs.org/en/latest/photutils/segmentation.html>`_.
Please see the Photutils `SourceProperties
<http://photutils.readthedocs.org/en/latest/api/photutils.segmentation.SourceProperties.html#photutils.segmentation.SourceProperties>`_
class for the list of the properties that are calculated for each
source.
