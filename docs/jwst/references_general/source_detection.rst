`~photutils.detection.DAOStarFinder` is an implementation of the `DAOFIND`_ algorithm
(`Stetson 1987, PASP 99, 191
<https://ui.adsabs.harvard.edu/abs/1987PASP...99..191S/abstract>`_).  It searches
images for local density maxima that have a peak amplitude greater
than a specified threshold (the threshold is applied to a convolved
image) and have a size and shape similar to a defined 2D Gaussian
kernel.  `DAOFIND`_ also provides an estimate of the object's
roundness and sharpness, whose lower and upper bounds can be
specified.

`~photutils.detection.IRAFStarFinder` is a Python implementation of the IRAF star finding algorithm,
which also calculates the objects' centroids, roundness, and sharpness.
However, `~photutils.detection.IRAFStarFinder` uses image moments
instead of 1-D Gaussian fits to projected light distributions like
`~photutils.detection.DAOStarFinder`.

`~photutils.segmentation.SourceFinder` implements an `image segmentation
<https://en.wikipedia.org/wiki/Image_segmentation>`_ algorithm, which is a
process of assigning a label to every pixel in an image such that
pixels with the same label are part of the same source.  The
segmentation procedure used is from
:ref:`Photutils source extraction <photutils:image_segmentation>`.
Detected sources must have a minimum number of connected pixels that
are each greater than a specified threshold value in an image.  The
threshold level is usually defined at some multiple of the background
standard deviation above the background.  The image can also be
filtered before thresholding to smooth the noise and maximize the
detectability of objects with a shape similar to the filter kernel.
Overlapping sources are detected as single sources.  Separating those
sources requires a deblending procedure, such as a multi-thresholding
technique used by `SExtractor
<https://www.astromatic.net/software/sextractor>`_.  Here we use the
`Photutils deblender
<https://photutils.readthedocs.io/en/stable/user_guide/segmentation.html#source-deblending>`_,
which is an algorithm that deblends sources using a combination of
multi-thresholding and `watershed segmentation
<https://en.wikipedia.org/wiki/Watershed_(image_processing)>`_.  In
order to deblend sources, they must be separated enough such that
there is a saddle between them.

.. warning::
    It has been shown (`STScI Technical Report JWST-STScI-008116, SM-12
    <https://www.stsci.edu/~goudfroo/NIRISSdoc/Centroid_Accuracies_Precisions_NIRISS_v2.pdf>`_)
    that for undersampled PSFs, e.g., for short-wavelength NIRISS
    imaging data, `~photutils.detection.DAOStarFinder` gives bad results no matter the input parameters
    due to its use of 1-D Gaussian fits.
    `~photutils.detection.IRAFStarFinder` or `~photutils.segmentation.SourceFinder` should be used instead.

.. _DAOFIND: https://ui.adsabs.harvard.edu/abs/1987PASP...99..191S/abstract
