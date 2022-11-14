Description
===========

:Classes: `jwst.resample.ResampleStep`, `jwst.resample.ResampleSpecStep`
:Alias: resample, resample_spec

This routine will resample each input 2D image based on the WCS and
distortion information, and will combine multiple resampled images
into a single undistorted product.  The distortion information should have
been incorporated into the image using the
:ref:`assign_wcs <assign_wcs_step>` step.

The ``resample`` step can take as input either:

  * a single 2D input image
  * an association table (in json format)

The defined parameters for the drizzle operation itself get
provided by the DRIZPARS reference file (from CRDS).  The exact values
used depends on the number of input images being combined and the filter
being used. Other information may be added as selection criteria later,
but for now, only basic information is used.

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

Context Image
-------------

In addition to image data, resample step also creates a "context image" stored
in the ``con`` attribute in the output data model or ``'CON'`` extension
of the FITS file. Each pixel in the context image is a bit field that encodes
information about which input image has contributed to the corresponding
pixel in the resampled data array. Context image uses 32 bit integers to encode
this information and hence it can keep track of only 32 input images.
First bit corresponds to the first input image, second bit corrsponds to the
second input image, and so on. If the number of input images is larger than 32,
then it is necessary to have multiple context images ("planes") to hold
information about all input images
with the first plane encoding which of the first 32 images contributed
to the output data pixel, second plane representing next 32 input images
(number 33-64), etc. For this reason, context array is a 3D array of the type
`numpy.int32` and shape ``(np, ny, nx)`` where ``nx`` and ``ny``
are dimensions of image's data. ``np`` is the number of "planes" equal to
``(number of input images - 1) // 32 + 1``. If a bit at position ``k`` in a
pixel with coordinates ``(p, y, x)`` is 0 then input image number
``32 * p + k`` (0-indexed) did not contribute to the output data pixel
with array coordinates ``(y, x)`` and if that bit is 1 then input image number
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
:py:func:`~jwst.resample.resample_utils.decode_context` function.

Spectroscopic Data
------------------

Use the ``resample_spec`` step for spectroscopic data.  The dispersion
direction is needed for this case, and this is obtained from the
DISPAXIS keyword.

References
----------

A full description of the drizzling algorithm can be found in
`Fruchter and Hook, PASP 2002 <https://doi.org/10.1086/338393>`_.
A description of the inverse variance map method can be found in
`Casertano et al., AJ 2000 <https://doi.org/10.1086/316851>`_, see Appendix A2.
A description of the drizzle parameters and other useful drizzle-related
resources can be found at `DrizzlePac Handbook <http://drizzlepac.stsci.edu>`_.
