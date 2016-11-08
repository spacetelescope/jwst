Reference File
==============
The reference file is a text file that uses JSON to hold the information
needed.

CRDS Selection Criteria
-----------------------
The file is selected based on the values of DETECTOR and FILTER (and
GRATING for NIRSpec).

Extract_1D Reference File Format
--------------------------------
All the information is specified in a list with key ``apertures``.  Each
element of this list is a dictionary, one for each aperture (e.g. a slit)
that is supported by the given reference file.  The particular dictionary
to use is found by matching the slit name in the science data with the
value of key ``id``.

The following keys are supported (but for IFU data, see below).
Key ``id`` is required for any element
of the ``apertures`` list that may be used; the value of ``id`` is compared
with the slit name (except for a full-frame input image) to select the
appropriate aperture.  Key ``dispaxis`` is similarly required.  Key
``region_type`` can be omitted, but if it is specified, its value must be
"target".  The source extraction region can be specified with ``ystart``,
``ystop``, etc., but a better alternative is to use ``src_coeff``.  If
background is to be subtracted, this should be specified by giving
``bkg_coeff``.  These are described in more detail below.  ``extract_width``
is not used if ``src_coeff`` is given.

* id: the slit name, e.g. "S200A1" (string)
* dispaxis: dispersion direction, 1 for X, 2 for Y (int)
* xstart: first pixel in X (int)
* xstop: last pixel in X (int)
* ystart: first pixel in Y (int)
* ystop: last pixel in Y (int)
* src_coeff: list of lists of float
* bkg_coeff: list of lists of float
* independent_var: "wavelength" or "pixel" (string)
* smoothing_width: width of boxcar for smoothing background regions along
  the dispersion direction (odd int)
* bkg_order: order of polynomial fit to background regions (int)
* extract_width: number of pixels in cross-dispersion direction (int)

For IFU data, these keys are used instead of most of the above:

* x_center: X pixel coordinate of the target (pixels, float)
* y_center: Y pixel coordinate of the target (pixels, float)
* extract_width: for a point source, this is the diameter of the circular
  extraction aperture; for an extended source, this is the width and height
  of the square extraction aperture (pixels, float)
* inner_bkg: (optional, and only for a point source) radius of the inner
  edge of the background annulus (pixels, float)
* outer_bkg: (optional, and only for a point source) radius of the outer
  edge of the background annulus (pixels, float)
* method: (optional, and only for a point source) one of "exact",
  "subpixel", or "center", the method used by photutils for computing the
  overlap between apertures and pixels (string, default is "exact")

Even if the extraction limits are specified by ``src_coeff`` (see below),
the limits in the dispersion direction can be specified by ``xstart`` and
``xstop`` (or ``ystart`` and ``ystop`` if ``dispaxis`` is 2).  These values
are zero-indexed pixel numbers within the subimage for a given "slit"; the
default is the size of the image in the dispersion direction.  The step
code may modify these values if, for example, the WCS information in the
input file has a more limited domain.

The source extraction region can be specified by giving ``src_coeff``,
coefficients for polynomial functions for the lower and upper limits of
the source extraction region.  Using this key will override the values
of ``ystart`` and ``ystop`` (if ``dispaxis`` is 1) or ``xstart`` and
``xstop`` (if ``dispaxis`` is 2), and ``extract_width``.  These polynomials
are functions of either wavelength (in microns) or pixel number (pixels in
the dispersion direction, with respect to the input 2-D slit image),
specified by the key ``independent_var`` (the default is "wavelength").
The values of these polynomial functions are pixel numbers in the
direction perpendicular to dispersion.  More than one source extraction
region may be specified, though this is not expected to be a typical case.

Background regions are specified by giving ``bkg_coeff``, coefficients for
polynomial functions for the lower and upper limits of one or more regions.
Background subtraction will be done only if ``bkg_coeff`` is given in the
reference file.  See below for an example.  See also ``bkg_order`` below.

The coefficients are specified as a list of an even number of lists (an
even number because both lower and upper limits must be specified).
The source extraction coefficients will normally be a list of just two
lists, the coefficients for the lower limit function and the coefficients
for the upper limit function of one extraction region.  The limits could
just be constant values, e.g. \[\[324.5\], \[335.5\]\].  Straight but tilted
lines are linear functions:

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
limits should be given as 324.5 and 335.5 respectively.

The order of a polynomial is specified implicitly to be one less than the
number of coefficients (this should not be confused with ``bkg_order``,
described below).  The number of coefficients must be at least one, and
there is no predefined upper limit.  The various polynomials (lower limits,
upper limits, possibly multiple regions) do not need to have the same
number of coefficients; each of the inner lists specifies a separate
polynomial.  However, the independent variable does need to be the same
for all polynomials for a given slit image (identified by key ``id``).

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
