Description
===========

This step creates a catalog of source photometry and morphologies.
Both aperture and isophotal (segment-based) photometry are calculated.
Source morphologies are based on 2D image moments within the source
segment.


Source Detection
^^^^^^^^^^^^^^^^

Sources are detected using `image segmentation
<https://en.wikipedia.org/wiki/Image_segmentation>`_, which is a
process of assigning a label to every pixel in an image such that
pixels with the same label are part of the same source.  The
segmentation procedure used is from `Photutils source extraction
<https://photutils.readthedocs.io/en/latest/segmentation.html>`_.
Detected sources must have a minimum number of connected pixels that
are each greater than a specified threshold value in an image.  The
threshold level is usually defined at some multiple of the background
standard deviation above the background.  The image can also be
filtered before thresholding to smooth the noise and maximize the
detectability of objects with a shape similar to the filter kernel.

Source Deblending
^^^^^^^^^^^^^^^^^

Overlapping sources are detected as single sources.  Separating those
sources requires a deblending procedure, such as a multi-thresholding
technique used by `SExtractor
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
photometry, centroids, and morphological properties.  The aperture
photometry is measured in three apertures, based on the input
encircled energy values.  The total aperture-corrected flux and
magnitudes are also calculated, based on the largest aperture.  Both
AB and Vega magnitudes are calculated.

For the isophotal photometry is based on `photutils segmentation
<https://photutils.readthedocs.org/en/latest/segmentation.html>`_.
The properties that are currently calculated for each source include
source centroids (both in pixel and sky coordinates), isophotal fluxes
(and errors), AB and Vega magnitudes (and errors), isophotal area,
semimajor and semiminor axis lengths, orientation of the major axis,
and sky coordinates at corners of the minimal bounding box enclosing
the source.

The output source catalog table is saved in `ECSV format
<https://docs.astropy.org/en/stable/io/ascii/write.html#ecsv-format>`_.
