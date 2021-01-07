Description
===========

Overview
--------
The ``extract_1d`` step extracts a 1D signal from a 2D or 3D dataset and
writes spectral data to an "x1d" product.  This works on all JWST spectroscopic
modes, including MIRI LRS (slit and slitless) and MRS, NIRCam WFSS and
TSGRISM, NIRISS WFSS and SOSS, and NIRSpec fixed-slit, IFU, and MOS.

An EXTRACT1D reference file is used for most modes to specify the location and
size of the target and background extraction apertures.
The EXTRACT1D reference file is not used for Wide-Field Slitless Spectroscopy data
(NIS_WFSS or NRC_WFSS). For these modes the extraction region is instead taken to be
the full size of the input 2D subarray or cutout for each source, or restricted to
the region within which the world coordinate system (WCS) is defined in each cutout.

For slit-like 2D input data, source and background extractions are done using
a rectangular aperture that covers one pixel in the dispersion direction and
uses a height in the cross-dispersion direction that is defined by parameters in
the EXTRACT1D reference file.
For 3D IFU data, on the other hand, the extraction options differ depending on
whether the target is a point or extended source.  For a point
source, the spectrum is extracted using circular aperture photometry,
optionally including background subtraction using a circular annulus.
For an extended source, rectangular aperture photometry is used, with
the entire image being extracted, and no background subtraction, regardless
of what was specified in the reference file or step arguments.
For either point or extended sources, the photometry makes use of the Astropy
affiliated package
`photutils <https://photutils.readthedocs.io/en/latest/>`_ to define an aperture
object and perform extraction.

For most spectral modes an aperture correction will be applied to the extracted
1D spectral data (unless otherwise selected by the user), in order to put the
results onto an infinite aperture scale.
This is done by creating interpolation functions based on the APCORR reference
file data and applying the interpolated aperture correction (a multiplicative
factor between 0 and 1) to the extracted, 1D spectral data (corrected data
include the "flux", "surf_bright", "error", and "sb_error" columns in the output
table).

Input
-----
Calibrated and potentially resampled 2D images or 3D cubes. The format should be a
CubeModel, SlitModel, IFUCubeModel, ImageModel, MultiSlitModel, or a ModelContainer.
For some JWST modes this is usually a resampled product, such as the "i2d" products
for MIRI LRS fixed-slit, NIRSpec fixed-slit, and NIRSpec MOS, or the "s3d" products
for MIRI MRS and NIRSpec IFU. For other modes that are not resampled (e.g. MIRI
LRS slitless, NIRISS SOSS, NIRSpec BrightObj, and NIRCam and NIRISS WFSS), this will
be a "cal" product.
For modes that have multiple slit instances (NIRSpec fixed-slit and MOS, WFSS),
The SCI extensions should have keyword SLTNAME to specify which slit was extracted,
though if there is only one slit (e.g. MIRI LRS and NIRISS SOSS), the slit name can
be taken from the EXTRACT1D reference file instead.

Normally the :ref:`photom <photom_step>` step should have been run before running
``extract_1d``.  If ``photom`` has not been run, a warning will be logged and the
output of ``extract_1d`` will be in units of count rate.  The ``photom`` step
converts data to units of either surface brightness (MegaJanskys per steradian) or,
for point sources observed with NIRSpec and NIRISS SOSS, units of flux density
(MegaJanskys).

Output
------
The output will be in ``MultiSpecModel`` format. For each input slit there will
be an output table extension with the name EXTRACT1D.  This extension will
have columns WAVELENGTH, FLUX, ERROR, SURF_BRIGHT, SB_ERROR, DQ,
BACKGROUND, BERROR and NPIXELS.
Some metadata will be written to the table header, mostly copied from the
input header.

The output WAVELENGTH data is copied from the wavelength array of the input 2D data,
if that attribute exists and was populated, otherwise it is calculated from the WCS.
FLUX is the flux density in Janskys; see keyword TUNIT2 if the data are
in a FITS BINTABLE.  ERROR is the error estimate for FLUX, and it has the
same units as FLUX.
SURF_BRIGHT is the surface brightness in MJy / sr, except that for point
sources observed with NIRSpec and NIRISS SOSS, SURF_BRIGHT will be set to
zero, because there's no way to express the extracted results from those modes
as a surface brightness. SB_ERROR is the error estimate for SURF_BRIGHT.
While it's expected that a user will make use of the FLUX column for
point-source data and the SURF_BRIGHT column for an extended source,
both columns are populated (except for NIRSpec and NIRISS SOSS point sources,
as mentioned above).
The ``extract_1d`` step collapses the input data from 2-D to 1-D by summing
one or more rows (or columns, depending on the dispersion direction).
A background may optionally be subtracted, but
there are also other options for background subtraction prior to ``extract_1d``.
For the case of input data in units of MJy / sr, the SURF_BRIGHT
and BACKGROUND columns are
populated by dividing the sum by the number of pixels (see the NPIXELS column,
described below) that were added together. The FLUX column is populated
by multiplying the sum by the solid angle of a pixel, and also multiplying
by 10^6 to convert from MJy to Jy.
For the case of input data in units of MJy (i.e. point sources,
NIRSpec or NIRISS SOSS), the SURF_BRIGHT column is set to zero, the
FLUX column is just multiplied by 10^6, and the BACKGROUND column is
divided by NPIXELS and by the solid angle of a pixel to convert to surface
brightness (MJy / sr).

NPIXELS is the number of pixels that were added together for the source
extraction region.  Note that this is not necessarily a constant, and
the value is not necessarily an integer (the data type is float).
BACKGROUND is the measured background, scaled to the extraction width used
for FLUX and SURF_BRIGHT.  BACKGROUND will be zero if background subtraction
is not requested.
ERROR, SB_ERROR, BERROR, and DQ are not populated with useful values yet.


.. _extract-1d-for-slits:

Extraction for 2D Slit Data
---------------------------
The operational details of the 1D extraction depend heavily on the parameter
values given in the :ref:`EXTRACT1D <extract1d_reffile>` reference file.
Here we describe their use within the ``extract_1d`` step.

Source Extraction Region
^^^^^^^^^^^^^^^^^^^^^^^^
As described in the documentation for the
:ref:`EXTRACT1D <extract1d_reffile>` reference file,
the characteristics of the source extraction region can be specified in one
of two different ways. 
The simplest approach is to use the ``xstart``, ``xstop``, ``ystart``,
``ystop``, and ``extract_width`` parameters.  Note that all of these values are
zero-indexed integers, and that the start and stop limits are inclusive.
If ``dispaxis=1``, the limits in the dispersion direction are ``xstart``
and ``xstop`` and the limits in the cross-dispersion direction are ``ystart``
and ``ystop``. If ``dispaxis=2``, the rolls are reversed.
If ``extract_width`` is also given, that takes priority over ``ystart`` and
``ystop`` (for ``dispaxis=1``) for the extraction width, but ``ystart`` and
``ystop`` will still be used to define the centering of the extraction region
in the cross-dispersion direction.
Any of these parameters will be modified internally by the step code if the
extraction region would extend outside the limits of the input image or outside
the domain specified by the WCS.

A more flexible way to specify the source extraction region is via the ``src_coeff``
parameter. ``src_coeff`` is specified as a list of lists of floating-point
polynomial coefficients that define the lower and upper
limits of the source extraction region as a function of dispersion. This allows,
for example, following a tilted or curved spectral trace or simply
following the variation in cross-dispersion FWHM as a function of wavelength.
If both ``src_coeff`` and ``ystart``/``ystop`` values are given, ``src_coeff``
takes precedence. The ``xstart`` and ``xstop`` values can still be used to
limit the range of the extraction in the dispersion direction. More details on
the specification and use of polynomial coefficients is given below.

Background Extraction Regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
One or more background extraction regions for a given aperture instance can
be specified using the ``bkg_coeff`` parameter in the EXTRACT1D reference file.
This is directly analogous to the use of ``src_coeff`` for specifiying source
extraction regions and functions in exactly the same way. More details on the
use of polynomial coefficients is given in the next section.
Background subtraction will be done if and only if ``bkg_coeff`` is given in
the EXTRACT1D reference file. The background is determined independently for
each column (or row, if dispersion is vertical), using pixel values from all
background regions within each column (or row).

Parameters related to background subtraction are ``smoothing_length``,
``bkg_fit``, and ``bkg_order``.

* If ``smoothing_length`` is specified, the 2D image data used to perform
  background extraction will be smoothed along the dispersion direction using
  a boxcar of width ``smoothing_length`` (in pixels). If not specified, no
  smoothing of the input 2D image data is performed.

* ``bkg_fit`` specifies the type of background computation to be performed
  within each column (or row). The default is "poly", in which case a
  polynomial of order ``bkg_order`` is fit to the background values within
  the column (or row). Alternatively, values of "mean" or "median" can be
  specified in order to compute the simple mean or median of the background
  values in each column (or row). Note that using "bkg_fit=mean" is
  mathematically equivalent to "bkg_fit=poly" with "bkg_order=0".

* If ``bkg_fit=poly`` is specified, ``bkg_order`` is used to indicate the
  polynomial order to be used. The default value is zero, i.e. a constant.

During source extraction, the background fit is evaluated at each pixel within the
source extraction region for that column (row), and the fitted values will
be subtracted (pixel by pixel) from the source count rate.

Source and Background Coefficient Lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The interpretation and use of polynomial coefficients to specify source and
background extraction regions via ``src_coeff`` and ``bkg_coeff`` is the same. 
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

Note: If the dispersion direction is vertical, replace "lower" with "left" and
"upper" with "right" in the above description.

Notice especially that ``src_coeff`` and ``bkg_coeff`` contain floating-point
values.  For interpreting fractions of a pixel, the convention used here
is that the pixel number at the center of a pixel is a whole number.  Thus,
if a lower or upper limit is a whole number, that limit splits the pixel
in two, so the weight for that pixel will be 0.5.  To include all the
pixels between 325 and 335 inclusive, for example, the lower and upper
limits would be given as 324.5 and 335.5 respectively.

The order of a polynomial is specified implicitly to be one less than the
number of coefficients. The number of coefficients for a lower or upper extraction
region limit must be at least one (i.e. zeroth-order polynomial). There is no
predefined upper limit on the number of coefficients (and hence polynomial order).
The various polynomials (lower limits, upper limits, possibly multiple regions) do
not need to have the same number of coefficients; each of the inner lists specifies
a separate polynomial. However, the independent variable (wavelength or pixel)
does need to be the same for all polynomials for a given slit.

Polynomials specified via ``src_coeff`` and ``bkg_coeff`` are functions of either wavelength
(in microns) or pixel number (pixels in the dispersion direction, with respect to
the input 2D slit image), which is specified by the parameter ``independent_var``.
The default is "pixel".  The values of these polynomial functions are pixel numbers in the
direction perpendicular to dispersion.

.. _extract-1d-for-ifu:

Extraction for 3D IFU Data
--------------------------
For IFU cube data, 1D extraction is contolled by a different set of EXTRACT1D
reference file parameters.
Note that for an extended source, anything specified in the reference file
or step arguments will be ignored; the entire image will be
extracted, and no background subtraction will be done.

For point sources a circular extraction aperture is used, along with an optional
circular annulus for background extraction and subtraction. The parameters
``x_center``, ``y_center``, and ``radius`` specify the location and size of the
source aperture, all of which are given in units of pixels.
The parameters ``inner_background`` and ``outer_background`` are used to
specify the limits of an annular background aperture. All radius-related
parameters can be floating-point values and appropriate partial-pixel scaling
will be used in the extraction process.

The region of overlap between an aperture and a pixel can be calculated by
one of three different methods, specified by the ``method`` parameter:  "exact"
(default), limited only by finite precision arithmetic; "center", the full value
in a pixel will be included if its center is within the aperture; or "subsample",
which means pixels will be subsampled N x N and the "center" option will be used
for each sub-pixel. When ``method`` is "subsample", the parameter ``subpixels``
is used to set the resampling value. The default value is 5.
