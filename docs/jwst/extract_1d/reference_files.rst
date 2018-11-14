Reference File
==============
The ``extract_1d`` step uses an EXTRACT1D reference file.

EXTRACT1D Reference File
------------------------
The EXTRACT1D reference file contains information needed to guide
the 1-D extraction process. It is a text file, with the information
in JSON format.

.. include:: extract1d_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for EXTRACT1D
++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in EXTRACT1D reference files,
because they are used as CRDS selectors
(see :ref:`extract1d_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
+++++++++++++++++++++
EXTRACT1D reference files are text files, with the information stored in
JSON format.
All the information is specified in a list with key ``apertures``.  Each
element of this list is a dictionary, one for each aperture (e.g. a slit)
that is supported by the given reference file.  The particular dictionary
to use is found by matching the slit name in the science data with the
value of key ``id``.  Key ``spectral_order`` is optional, but if it is
present, it must match the expected spectral order number.

The following keys are supported for non-IFU data (see below for IFU keys).
Key ``id`` is the primary criterion for selecting which element of
the ``apertures`` list to use.  The slit name (except for a full-frame
input image) is compared with the values of ``id`` in the ``apertures``
list to select the appropriate aperture.
In order to allow the possibility of multiple
spectral orders for the same slit name, there may be more than one element
of ``apertures`` with the same value for key ``id``.  These should then be
distinguished by using the secondary selection criterion ``spectral_order``.
In this case, the various spectral orders would likely have different
extraction locations within the image, so different elements of ``apertures``
are needed in order to specify those locations.
If key ``dispaxis`` is specified, its value will be used to set the
dispersion direction within the image.  If ``dispaxis`` is
not specified, the dispersion direction will be taken to be the axis
along which the wavelengths change more rapidly.
Key ``region_type`` can be omitted, but if it is specified, its value must
be "target".  The source extraction region can be specified with ``ystart``,
``ystop``, etc., but a more flexible alternative is to use ``src_coeff``.
If background is to be subtracted, this should be specified by giving
``bkg_coeff``.  These are described in more detail below.

* id: the slit name, e.g. "S200A1" (string)
* spectral_order: the spectral order number (optional); this can be either
  positive or negative, but it should not be zero (int)
* dispaxis: dispersion direction, 1 for X, 2 for Y (int)
* xstart: first pixel in the horizontal direction, X (int)
* xstop: last pixel in the horizontal direction, X (int)
* ystart: first pixel in the vertical direction, Y (int)
* ystop: last pixel in the vertical direction, Y (int)
* src_coeff: this takes priority for specifying the source extraction region
  (list of lists of float)
* bkg_coeff: for specifying background subraction regions
  (list of lists of float)
* independent_var: "wavelength" or "pixel" (string)
* smoothing_length: width of boxcar for smoothing background regions along
  the dispersion direction (odd int)
* bkg_order: order of polynomial fit to background regions (int)
* extract_width: number of pixels in cross-dispersion direction (int)

If ``src_coeff`` is given, those coefficients take priority for specifying
the source extraction region in the cross-dispersion direction.  ``xstart``
and ``xstop`` (or ``ystart`` and ``ystop`` if ``dispaxis`` is 2) will
still be used for the limits in the dispersion direction.  Background
subtraction will be done if and only if ``bkg_coeff`` is given.  See below
for further details.

For IFU cube data, the following keys are used instead of those above:

* id: the slit name, but this can be "ANY" (string)
* x_center: X pixel coordinate of the target (pixels, float, the default
  is the center of the image along the X axis)
* y_center: Y pixel coordinate of the target (pixels, float, the default
  is the center of the image along the Y axis)
* radius: (only used for a point source) the radius of the circular
  extraction aperture (pixels, float, default is one quarter of the smaller
  of the image axis lengths)
* subtract_background: (only used for a point source) if true, subtract a
  background determined from an annulus with inner and outer radii given
  by ``inner_bkg`` and ``outer_bkg`` (boolean)
* inner_bkg: (only for a point source) radius of the inner edge of the
  background annulus (pixels, float, default = ``radius``)
* outer_bkg: (only for a point source) radius of the outer edge of the
  background annulus (pixels, float, default = ``inner_bkg * sqrt(2)``)
* width: (only for an extended source) the width of the rectangular
  extraction region; if ``theta = 0``, the width side is along the X axis
  (pixels, float, default is half of the smaller image axis length)
* height: (only for an extended source) the height of the rectangular
  extraction region; if ``theta = 0``, the height side is along the Y axis
  (pixels, float, default is half of the smaller image axis length)
* angle: (only for an extended source) the counterclockwise rotation angle of
  the ``width`` side from the positive X axis (degrees)
* method: one of "exact", "subpixel", or "center", the method
  used by photutils for computing the overlap between apertures and pixels
  (string, default is "exact")
* subpixels: if ``method`` is "subpixel", pixels will be resampled by this
  factor in each dimension (int, the default is 5)

The rest of this description pertains to the parameters for non-IFU data.

If ``src_coeff`` is not given, the extraction limits can be specified by
``xstart``, ``xstop``, ``ystart``, ``ystop``, and ``extract_width``.  Note
that all of these values are integers, and that the start and stop limits
are inclusive.
If ``dispaxis``
is 1, the zero-indexed limits in the dispersion direction are ``xstart``
and ``xstop``; if ``dispaxis`` is 2, the dispersion limits are ``ystart``
and ``ystop``.  (The dispersion limits can be given even if ``src_coeff``
has been used for defining the cross-dispersion limits.)  The limits in
the cross-dispersion direction can be given by ``ystart`` and ``ystop``
(or ``xstart`` and ``xstop`` if ``dispaxis`` is 2).  If ``extract_width``
is also given, that takes priority over ``ystart`` to ``ystop`` (for
``dispaxis`` = 1) for the extraction width, but ``ystart`` and ``ystop``
(for ``dispaxis`` = 1) will still be used to define the middle in the
cross-dispersion direction.  Any of these parameters can be modified
by the step code if the extraction region would extend outside the input
image, or outside the domain specified by the WCS.

The source extraction region can be specified more precisely by giving
``src_coeff``, coefficients for polynomial functions for the lower and
upper limits of the source extraction region.  As described in the previous
paragraph, using this key will override the values
of ``ystart`` and ``ystop`` (if ``dispaxis`` is 1) or ``xstart`` and
``xstop`` (if ``dispaxis`` is 2), and ``extract_width``.  These polynomials
are functions of either wavelength (in microns) or pixel number (pixels in
the dispersion direction, with respect to the input 2-D slit image),
specified by the key ``independent_var``.  The default is "pixel".
The values of these polynomial functions are pixel numbers in the
direction perpendicular to dispersion.  More than one source extraction
region may be specified, though this is not expected to be a typical case.

Background regions are specified by giving ``bkg_coeff``, coefficients for
polynomial functions for the lower and upper limits of one or more regions.
Background subtraction will be done only if ``bkg_coeff`` is given in the
reference file.  See below for an example.  See also ``bkg_order`` below.

The coefficients are specified as a list of an even number of lists (an
even number because both the lower and upper limits of each extraction region
must be specified).  The source extraction coefficients will normally be
a list of just two lists, the coefficients for the lower limit function
and the coefficients for the upper limit function of one extraction
region.  The limits could just be constant values,
e.g. \[\[324.5\], \[335.5\]\].  Straight but tilted lines are linear functions:

\[\[324.5, 0.0137\], \[335.5, 0.0137\]\]

Multiple regions may be specified for either the source or background, or
both.  It will be common to specify more than one background region.  Here
is an example for specifying two background regions:

\[\[315.2, 0.0135\], \[320.7, 0.0135\], \[341.1, 0.0139\], \[346.8, 0.0139\]\]

This is interpreted as follows:

* \[315.2, 0.0135\]: lower limit for first background region
* \[320.7, 0.0135\]: upper limit for first background region
* \[341.1, 0.0139\]: lower limit for second background region
* \[346.8, 0.0139\]: upper limit for second background region

If the dispersion direction is vertical, replace "lower" with "left" and
"upper" with "right" in the above description.

Note especially that ``src_coeff`` and ``bkg_coeff`` contain floating-point
values.  For interpreting fractions of a pixel, the convention used here
is that the pixel number at the center of a pixel is a whole number.  Thus,
if a lower or upper limit is a whole number, that limit splits the pixel
in two, so the weight for that pixel will be 0.5.  To include all the
pixels between 325 and 335 inclusive, for example, the lower and upper
limits would be given as 324.5 and 335.5 respectively.

The order of a polynomial is specified implicitly to be one less than the
number of coefficients (this should not be confused with ``bkg_order``,
described below).  The number of coefficients must be at least one, and
there is no predefined upper limit.  The various polynomials (lower limits,
upper limits, possibly multiple regions) do not need to have the same
number of coefficients; each of the inner lists specifies a separate
polynomial.  However, the independent variable (wavelength or pixel)
does need to be the same for all polynomials for a given slit image
(identified by key ``id``).

The background is determined independently for each column (or row, if
``dispaxis`` is 2) of the spectrum.  The ``smoothing_length`` parameter
is the width of a boxcar for smoothing the background in the dispersion
direction.  If this is not specified, either in the reference file, the
config file, or on the command line, no smoothing will be done along the
dispersion direction.  Following background smoothing (if any), for each
column (row), a polynomial of order ``bkg_order`` will be fit to the pixel
values in that column (row) in all the background regions.  If not
specified, a value of 0 will be used, i.e. a constant function, the mean
value.  The polynomial will then be evaluated at each pixel within the
source extraction region for that column (row), and the fitted values will
be subtracted (pixel by pixel) from the source count rate.

Example EXTRACT1D Reference File
--------------------------------
The following JSON was taken as an example from reference file
jwst_niriss_extract1d_0003.json::
  {
      "REFTYPE": "EXTRACT1D",
      "INSTRUME": "NIRISS",
      "TELESCOP": "JWST",
      "DETECTOR": "NIS",
      "EXP_TYPE": "NIS_SOSS",
      "PEDIGREE": "GROUND",
      "DESCRIP": "NIRISS SOSS extraction params for ground testing",
      "AUTHOR": "M.Wolfe, H.Bushouse",
      "HISTORY": "This reference file is for used in Build 7.1 of the JWST Calibraton pipeline. The regions are rectagular and do not follow the trace.",
      "USEAFTER": "2015-11-01T00:00:00",
      "apertures": [
        {
        "id": "FULL",
        "region_type": "target",
        "bkg_coeff": [[2014.5],[2043.5]],
        "xstart":    4,
        "xstop":  2044,
        "ystart": 1792,
        "ystop":  1972,
        "dispaxis": 1,
        "extract_width": 181
        },

        {
        "id": "SUBSTRIP256",
        "region_type": "target",
        "bkg_coeff": [[221.5],[251.5]],
        "xstart":   4,
        "xstop": 2044,
        "ystart": 20,
        "ystop":  220,
        "dispaxis": 1,
        "extract_width": 201
        },

        {
        "id": "SUBSTRIP96",
        "region_type": "target",
        "bkg_coeff": [[1.5],[8.5],[92.5],[94.5]],
        "xstart":   4,
        "xstop": 2044,
        "ystart":  10,
        "ystop":   92,
        "dispaxis": 1,
        "extract_width": 83
        }]
  }

