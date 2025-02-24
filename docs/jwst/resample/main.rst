Description
===========

:Classes: `jwst.resample.ResampleStep`
:Alias: resample

This routine will resample each input 2D image based on the WCS and
distortion information, and will combine multiple resampled images
into a single undistorted product.  The distortion information should have
been incorporated into the image using the
:ref:`assign_wcs <assign_wcs_step>` step.

The ``resample`` step can take as input either:

#. a single 2D input image
#. an association table (in json format)

The parameters for the drizzle operation are provided via ``resample``
step parameters, which can be overridden by a step parameter reference
file from CRDS.

The output product gets defined using the WCS information of all inputs,
even if it is just a single input image. The output WCS defines a
field-of-view that encompasses the undistorted footprints on the sky
of all the input images with the same orientation and plate scale
as the first listed input image.

This step uses the interface to the C-based cdriz routine to do the
resampling via the drizzle method.  The input-to-output pixel
mapping is determined via a mapping function derived from the
WCS of each input image and the WCS of the defined output product.
This mapping function gets passed to cdriz to drive the actual
drizzling to create the output product.


Error Propagation
-----------------

The error associated with each resampled pixel can in principle be derived
from the variance components associated with each input pixel, weighted by
the square of the input user weights and the square of the overlap between
the input and output pixels. In practice, the cdriz routine does not currently
support propagating variance data alongside science images, so the output
error cannot be precisely calculated.

To approximate the error on a resampled pixel, the variance arrays associated
with each input image are resampled individually, then combined with a weighted
sum.  The process is:

#. For each input image, take the square root of each of the read noise variance
   array to make an error image.

#. Drizzle the read noise error image onto the output WCS, with drizzle
   parameters matching those used for the science data.

#. Square the resampled read noise to make a variance array.

   a. If the resampling `weight_type` is an inverse variance map (`ivm`), weight
      the resampled variance by the square of its own inverse.

   #. If the `weight_type` is the exposure time (`exptime`), weight the
      resampled variance by the square of the exposure time for the image.

#. Add the weighted, resampled read noise variance to a running sum across all
   images.  Add the weights (unsquared) to a separate running sum across
   all images.

#. Perform the same steps for the Poisson noise variance and the flat variance.
   For these components, the weight for the sum is either the resampled read
   noise variance or else the exposure time.

#. For each variance component (read noise, Poisson, and flat), divide the
   weighted variance sum by the total weight, squared.

After each variance component is resampled and summed, the final error
array is computed as the square root of the sum of the three independent
variance components.  This error image is stored in the ``err`` attribute
in the output data model or the ``'ERR'`` extension of the FITS file.

It is expected that the output errors computed in this way will
generally overestimate the true error on the resampled data.  The magnitude
of the overestimation depends on the details of the pixel weights
and error images.  Note, however, that drizzling error images produces
a much better estimate of the output error than directly drizzling
the variance images, since the kernel overlap weights do not need to be
squared for combining error values.


Context Image
-------------

In addition to image data, resample step also creates a "context image" stored
in the ``con`` attribute in the output data model or ``'CON'`` extension
of the FITS file. Each pixel in the context image is a bit field that encodes
information about which input image has contributed to the corresponding
pixel in the resampled data array. Context image uses 32 bit integers to encode
this information and hence it can keep track of only 32 input images.
First bit corresponds to the first input image, second bit corresponds to the
second input image, and so on. If the number of input images is larger than 32,
then it is necessary to have multiple context images ("planes") to hold
information about all input images,
with the first plane encoding which of the first 32 images contributed
to the output data pixel, the second plane representing next 32 input images
(number 33-64), etc. For this reason, context array is a 3D array of the type
`numpy.int32` and shape ``(np, ny, nx)`` where ``nx`` and ``ny``
are dimensions of the image data. ``np`` is the number of "planes" computed as
``(number of input images - 1) // 32 + 1``. If a bit at position ``k`` in a
pixel with coordinates ``(p, y, x)`` is 0, then input image number
``32 * p + k`` (0-indexed) did not contribute to the output data pixel
with array coordinates ``(y, x)`` and if that bit is 1, then input image number
``32 * p + k`` did contribute to the pixel ``(y, x)`` in the resampled image.

As an example, let's assume we have 8 input images. Then, when ``'CON'`` pixel
values are displayed using binary representation (and decimal in parenthesis),
one could see values like this::

    00000001 (1) - only first input image contributed to this output pixel;
    00000010 (2) - 2nd input image contributed;
    00000100 (4) - 3rd input image contributed;
    10000000 (128) - 8th input image contributed;
    10000100 (132=128+4) - 3rd and 8th input images contributed;
    11001101 (205=1+4+8+64+128) - input images 1, 3, 4, 7, 8 have contributed
    to this output pixel.

In order to test if a specific input image contributed to an output pixel,
one needs to use bitwise operations. Using the example above, to test whether
input images number 4 and 5 have contributed to the output pixel whose
corresponding ``'CON'`` value is 205 (11001101 in binary form) we can do
the following:

>>> bool(205 & (1 << (5 - 1)))  # (205 & 16) = 0 (== 0 => False): did NOT contribute
False
>>> bool(205 & (1 << (4 - 1)))  # (205 & 8) = 8 (!= 0 => True): did contribute
True

In general, to get a list of all input images that have contributed to an
output resampled pixel with image coordinates ``(x, y)``, and given a
context array ``con``, one can do something like this:

.. doctest-skip::

    >>> import numpy as np
    >>> np.flatnonzero([v & (1 << k) for v in con[:, y, x] for k in range(32)])

For convenience, this functionality was implemented in the
:py:func:`~drizzle.utils.decode_context` function.


References
----------

A full description of the drizzling algorithm can be found in
`Fruchter and Hook, PASP 2002 <https://doi.org/10.1086/338393>`_.
A description of the inverse variance map method can be found in
`Casertano et al., AJ 2000 <https://doi.org/10.1086/316851>`_, see Appendix A2.
A description of the drizzle parameters and other useful drizzle-related
resources can be found at `DrizzlePac Handbook <http://drizzlepac.stsci.edu>`_.
