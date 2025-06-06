Description
===========

:Class: `jwst.extract_1d.Extract1dStep`
:Alias: extract_1d

Overview
--------
The ``extract_1d`` step extracts a 1D signal from a 2D or 3D dataset and
writes spectral data to an "x1d" product (or "x1dints" for time series data).
This step works on all JWST spectroscopic modes, including MIRI LRS (slit and slitless)
and MRS, NIRCam WFSS and TSGRISM, NIRISS WFSS and SOSS, and NIRSpec fixed-slit, IFU, and MOS.

An EXTRACT1D reference file is used for most modes to specify the location and
size of the target and background extraction apertures.
The EXTRACT1D reference file is not used for Wide-Field Slitless Spectroscopy data
(NIS_WFSS or NRC_WFSS). For these modes the extraction region is instead taken to be
the full size of the input 2D subarray or cutout for each source, or restricted to
the region within which the world coordinate system (WCS) is defined in each cutout.

For slit-like 2D input data, source and background extractions are, by default, done
using a rectangular aperture that covers one pixel in the dispersion direction and
uses a height in the cross-dispersion direction that is defined by parameters in
the EXTRACT1D reference file.  Optionally, for point sources, a PSF-based optimal
extraction may be performed, using a model of the spectral PSF to fit the total flux
at each dispersion element.

For 3D IFU data, the extraction options differ depending on
whether the target is a point or extended source.  For a point
source, the spectrum is extracted using circular aperture photometry,
optionally including background subtraction using a circular annulus.
For an extended source, rectangular aperture photometry is used, with
the entire image being extracted, and no background subtraction, regardless
of what was specified in the reference file or step arguments.
For both point or extended sources, photometric measurements make use of
the Astropy affiliated package
`photutils <https://photutils.readthedocs.io/en/latest/>`_ to define an aperture
object and perform extraction.  For 3D NIRSpec fixed slit rateints data, the
``extract_1d`` step will be skipped as 3D input for the mode is not supported.


For most spectral modes, an aperture correction will be applied to the extracted
1D spectral data (unless otherwise selected by the user), in order to put the
results onto an infinite aperture scale.
This is done by creating interpolation functions based on the APCORR reference
file data and applying the interpolated aperture correction (a multiplicative
factor between 0 and 1) to the extracted, 1D spectral data (corrected data
include the "flux", "surf_bright", "flux_error", "sb_error", and all flux and
surface brightness variance columns in the output table).  For optimal extractions,
aperture correction is not performed, since it is assumed the total flux is
modeled by the PSF.

Input
-----
The input data are calibrated and potentially resampled 2D images or 3D cubes.
The format should be a
CubeModel, SlitModel, IFUCubeModel, ImageModel, MultiSlitModel, or a ModelContainer.
For some JWST modes this is usually a resampled product, such as the "s2d" products
for MIRI LRS fixed-slit, NIRSpec fixed-slit, and NIRSpec MOS, or the "s3d" products
for MIRI MRS and NIRSpec IFU. For other modes that are not resampled (e.g. MIRI
LRS slitless, NIRISS SOSS, NIRSpec BOTS, and NIRCam and NIRISS WFSS), this will
be a "cal" or "calints" product.
For modes that have multiple slit instances (NIRSpec fixed-slit and MOS, WFSS),
the SCI extensions should have the keyword SLTNAME to specify which slit was extracted.

Normally the :ref:`photom <photom_step>` step should be applied before running
``extract_1d``.  If ``photom`` has not been run, a warning will be logged and the
output of ``extract_1d`` will be in units of count rate.  The ``photom`` step
converts data to units of either surface brightness (megajanskys per steradian) or,
for point sources observed with NIRSpec and NIRISS SOSS, units of flux density
(megajanskys).

Output
------

Output data for this step are in the form of :ref:`x1d`.

Data structure
^^^^^^^^^^^^^^

The output for most modes will be in ``MultiSpecModel`` format. This datamodel collects
multiple spectra in a list, stored in the ``spec`` attribute.  The data for each spectrum
is stored in a table under the ``spec_table`` attribute for the spectrum.  The spectral models
in the ``spec`` attribute also hold related metadata, such as the slit name (``name``) or the
source ID (``source_id``).

In the output file, the spectral data is stored as a table extension with the name EXTRACT1D.
This extension will have columns WAVELENGTH, FLUX, FLUX_ERROR, FLUX_VAR_POISSON, FLUX_VAR_RNOISE,
FLUX_VAR_FLAT, SURF_BRIGHT, SB_ERROR, SB_VAR_POISSON, SB_VAR_RNOISE,
SB_VAR_FLAT, DQ, BACKGROUND, BKGD_ERROR, BKGD_VAR_POISSON, BKGD_VAR_RNOISE,
BKGD_VAR_FLAT and NPIXELS.

For example, to access the slit name, wavelength, and flux from each spectrum in a model:

.. doctest-skip::

  >>> from stdatamodels.jwst import datamodels
  >>> multi_spec = datamodels.open('multi_spec_x1d.fits')
  >>> for spectrum in multi_spec.spec:
  >>>     slit_name = spectrum.name
  >>>     wave = spectrum.spec_table["WAVELENGTH"]
  >>>     flux = spectrum.spec_table["FLUX"]


In the case of MIRI MRS data, the output is a ``MRSMultiSpecModel``. This model has the
same structure as a ``MultiSpecModel``, except that there are three additional
columns in the output table:  RF_FLUX, RF_SURF_BRIGHT, and RF_BACKGROUND.
For more details on the MIRI MRS extracted data see :ref:`MIRI-MRS-1D-residual-fringe`.


For NIRCam and NIRISS WFSS data, hundreds to thousands of spectra from different sources
may be extracted. For those modes, the output is a ``WFSSMultiSpecModel``.
The data in this model is stored in the ``spec`` attribute, such that one spectral table
is created for each exposure for each spectral order in the input data.
Each extension of the output FITS file thus represents one exposure/spectral order combination.
The extension metadata contains a unique exposure ID (FITS keyword EXPGRPID) for each extension,
which combines exposure grouping metadata, the exposure number, and the spectral order.
The spectral tables for this model contain the same columns as the ``MultiSpecModel``, but
each row in the table contains the full spectrum for a single source and order. The spectral columns
are 2D: each row is a 1D vector containing all data points for the spectrum. In addition, the
spectral tables for this model have extra 1D columns to contain the metadata for the spectrum in each row.
These metadata fields include:
SOURCE_ID, N_ALONGDISP, SOURCE_TYPE, SOURCE_XPOS, SOURCE_YPOS, SOURCE_RA, SOURCE_DEC, 
EXTRACT2D_XSTART, EXTRACT2D_YSTART.

Note that the vector columns have the same length for all the sources in the table, meaning that
the number of elements in the table rows is set by the spectrum with the most data points.
The other spectra are NaN-padded to match the longest spectrum,
and the number of valid data points for each spectrum is recorded
in the N_ALONGDISP column.

For example, to access the wavelength and flux for a specific source ID (say, 1200) and
integration (the first) in a WFSSMultiSpecModel:

.. doctest-skip::

  >>> from stdatamodels.jwst import datamodels
  >>> model = datamodels.open('multi_wfss_x1d.fits')
  >>> first_table = model.spec[0].spec_table
  >>> id_want = 1200
  >>> row_want = first_table[first_table["SOURCE_ID"] == id_want][0]
  >>> nelem = row_want["N_ALONGDISP"]
  >>> wave, flux = row_want["WAVELENGTH"][:nelem], row_want["FLUX"][:nelem]


For time series observations (TSO) with spectra extracted from multiple integrations,
the output is a ``TSOMultiSpecModel``.  The spectral tables for this model have
the same columns as the ``MultiSpecModel``, but each row in the table contains the full
spectrum for a single integration.  The spectral columns are 2D: each row is a 1D
vector containing all data points for the spectrum in that integration.  The spectral tables
for this model have extra 1D columns to contain the metadata for the spectrum in each row.
These metadata fields include the segment number, integration number, and various time tags,
as follows:
SEGMENT, INT_NUM, START_TIME_MJD, MID_TIME_MJD, END_TIME_MJD, START_TDB, MID_TDB, and END_TDB.

For example, to access the slit name, integration number, wavelength, and flux from
each spectrum in a TSO model:

.. doctest-skip::

  >>> from stdatamodels.jwst import datamodels
  >>> multi_int_spec = datamodels.open('multi_spec_x1dints.fits')
  >>> for spectrum in multi_int_spec.spec:
  >>>     slit_name = spectrum.name
  >>>     integrations = spectrum.spec_table["INT_NUM"]
  >>>     for i, int_num in enumerate(integrations):
  >>>         wave = spectrum.spec_table["WAVELENGTH"][i]
  >>>         flux = spectrum.spec_table["FLUX"][i]


Data sources
^^^^^^^^^^^^

For all modes, some metadata for the slit and source will be written to the header for
the table extension, mostly copied from the input SCI extension headers.

For slit-like modes, the extraction region is
recorded in the metadata of the table header as EXTRXSTR (x start of extraction),
EXTRXSTP (x end of extraction), EXTRYSTR (y start of extraction), and
EXTRYSTP (y end of extraction).  For MIRI and NIRSpec IFU data, the center of
the extraction region is recorded in the metadata EXTR_X (x center of extraction region)
and EXTR_Y (y center of extraction region). The NIRISS SOSS algorithm is a specialized extraction
algorithm that does not use fixed limits; therefore, no extraction limits are provided for this mode.
Note that the pipeline takes input start/stop values from the reference files to be
zero-indexed positions, but all extraction values are recorded in the headers as one-indexed
values, following FITS header conventions.

The output WAVELENGTH data is copied from the wavelength array of the input 2D data,
if that attribute exists and was populated. Otherwise, it is calculated from the WCS.
FLUX is the summed flux density in janskys (see keyword TUNIT2 in the FITS table header).
FLUX_ERROR is the error estimate for FLUX; it has the
same units as FLUX. The error is calculated as the square root of the sum of the
three variance arrays: Poisson, read noise (RNOISE), and flat field (FLAT).
SURF_BRIGHT is the surface brightness in MJy / sr, except that for point
sources observed with NIRSpec and NIRISS SOSS, or optimal extractions, SURF_BRIGHT will be set to
zero, because there is no way to express the extracted results from those modes
as a surface brightness. SB_ERROR is the error estimate for SURF_BRIGHT, calculated
in the same fashion as FLUX_ERROR but using the SB_VAR arrays. While it's expected
that a user will make use of the FLUX column for point-source data and the
SURF_BRIGHT column for an extended source, both columns are populated
(except as mentioned above).

The ``extract_1d`` step collapses the input data from 2-D to 1-D by summing
one or more rows (or columns, depending on the dispersion direction).
A residual background may optionally be subtracted, in addition to any
background subtraction performed prior to ``extract_1d``.
For the case of input data in units of MJy / sr, the SURF_BRIGHT
and BACKGROUND columns are populated by dividing the sum by the number of pixels
(see the NPIXELS column, described below) summed over during extraction.
The FLUX column is populated by multiplying the sum by the solid angle of a pixel,
and also multiplying by 10^6 to convert from MJy to Jy.
For the case of input data in units of MJy (i.e. point sources,
NIRSpec or NIRISS SOSS), the SURF_BRIGHT column is set to zero, the
FLUX column is just multiplied by 10^6, and the BACKGROUND column is
divided by NPIXELS and by the solid angle of a pixel to convert to surface
brightness (MJy / sr).

NPIXELS is the number of pixels that were added together for the source
extraction region.  Note that this is not necessarily a constant, since some
pixels might be excluded for some wavelengths and included for others, and
the value is not necessarily an integer, since partial pixels may have been
included in the extraction aperture.

BACKGROUND is the measured background, scaled to the extraction width used
for FLUX and SURF_BRIGHT.  BACKGROUND will be zero if background subtraction
is not requested. BKGD_ERROR is calculated as the square root of the sum of the
BKGD_VAR arrays.

The DQ array is set to DO_NOT_USE for pixels with NaN flux values and zero
otherwise.


.. _extract-1d-for-slits:

Box Extraction for 2D Slit Data
-------------------------------
For standard box extractions, the operational details depend heavily on the parameter
values given in the :ref:`EXTRACT1D <extract1d_reffile>` reference file.
Here we describe their use within the ``extract_1d`` step.

Source Extraction Region
^^^^^^^^^^^^^^^^^^^^^^^^
As described in the documentation for the
:ref:`EXTRACT1D <extract1d_reffile>` reference file,
the characteristics of the source extraction region can be specified in one
of two different ways.

The simplest approach is to use the `xstart`, `xstop`, `ystart`,
`ystop`, and `extract_width` parameters.  Note that all of these values are
zero-indexed floating point values, the start and stop limits are inclusive, and
the values are in the frame of the image being operated on (which could be a cutout
of a larger original image).
If `dispaxis=1`, the limits in the dispersion direction are `xstart`
and `xstop` and the limits in the cross-dispersion direction are `ystart`
and `ystop`. If `dispaxis=2`, the roles are reversed.

If `extract_width` is also given, the start and stop values are used to define
the center of the extraction region in the cross-dispersion direction, but the
width of the aperture is set by the `extract_width` value.

For some instruments and modes, the extraction region may be adjusted
to account for the expected location of the source with the `use_source_posn` 
option. This option is available for NIRSpec MOS, fixed-slit, and BOTS data, 
as well as MIRI LRS fixed-slit.
If `use_source_posn` is set to None via the reference file or input parameters,
it is turned on by default for all point sources in these modes.
To turn it on for extended sources, set `use_source_posn` to True.
To turn it off for any mode, set `use_source_posn` to False.
If source position option is enabled, the planned location for the source and its 
trace are calculated internally via header metadata recording the source position 
and the spectral WCS transforms.  The source location will be used to offset the 
extraction start and stop values in the cross-dispersion direction.
If `extract_width` is provided, the source extraction region will be centered 
on the calculated trace with a width set by the `extract_width` value.  
For resampled, "s2d", products this will effectively be the rectangular 
extraction region offset in the cross-dispersion direction.  For 
"cal" or "calints" products that have not been resampled, the extraction region 
will be curved to follow the calculated trace.
If no `extract_width` has been provided, the shifted extraction start and 
stop values will be used.

A more flexible way to specify the source extraction region is via the `src_coeff`
parameter. `src_coeff` is specified as a list of lists of floating-point
polynomial coefficients that define the lower and upper
limits of the source extraction region as a function of dispersion. This allows,
for example, following a tilted or curved spectral trace or simply
following the variation in cross-dispersion FWHM as a function of wavelength.
If both `src_coeff` and cross-dispersion start/stop values are given, `src_coeff`
takes precedence. The start/stop values can still be used to
limit the range of the extraction in the dispersion direction. More details on
the specification and use of polynomial coefficients is given below.

Note that if source position correction is enabled, the position offset is applied to
any supplied `src_coeff` values, as well as the cross-dispersion start/stop values.
To ensure the provided `src_coeff` values are used as-is, set `use_source_posn`
to False.


Background Extraction Regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
One or more background extraction regions for a given aperture instance can
be specified using the `bkg_coeff` parameter in the EXTRACT1D reference file.
This is directly analogous to the use of `src_coeff` for specifying source
extraction regions and functions in exactly the same way. More details on the
use of polynomial coefficients is given in the next section.

By default, background subtraction will be done if `bkg_coeff` is set in
the EXTRACT1D reference file. To turn it off without modifying the reference
file, set `subtract_background` to False in the input step parameters.

The background values are determined independently for
each column (or row, if dispersion is vertical), using pixel values from all
background regions within each column (or row).
Parameters related to background fitting are `smoothing_length`,
`bkg_fit`, and `bkg_order`:

#. If `smoothing_length` is specified, the 2D image data used to perform
   background extraction will be smoothed along the dispersion direction using
   a boxcar of width `smoothing_length` (in pixels). If not specified, no
   smoothing of the input 2D image data is performed.

#. `bkg_fit` specifies the type of fit to the background data, to be performed
   within each column (or row). The default value is None; if not set by
   the user, the step will search the reference file for a value. If no value
   is found, `bkg_fit` will be set to "poly". The "poly" mode fits a
   polynomial of order `bkg_order` to the background values within
   the column (or row). Alternatively, values of "mean" or "median" can be
   specified in order to compute the simple mean or median of the background
   values in each column (or row). Note that using `bkg_fit=mean` is
   mathematically equivalent to `bkg_fit=poly` with `bkg_order=0`.

#. If `bkg_fit=poly` is specified, `bkg_order` is used to indicate the
   polynomial order to be used. The default value is zero, i.e. a constant.

During source extraction, the background fit is evaluated at each pixel within the
source extraction region for that column/row, and the fitted values will
be subtracted (pixel by pixel) from the source count rate, prior to summing
over the aperture.

If source position correction is enabled, the calculated position offset is applied to
any supplied `bkg_coeff` values, as well as the source aperture limit values.
To ensure the provided `bkg_coeff` values are used as-is, set `use_source_posn`
to False.

Source and Background Coefficient Lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The interpretation and use of polynomial coefficients to specify source and
background extraction regions is the same for both source coefficients (`src_coeff`)
and background coefficients (`bkg_coeff`).

Polynomials specified via `src_coeff` and `bkg_coeff` are functions of either wavelength
(in microns) or pixel number (pixels in the dispersion direction, with respect to
the input 2D slit image), which is specified by the parameter `independent_var`.
The default is "pixel"; the alternative is "wavelength".  The dependent values of these
polynomial functions are always pixel numbers (zero-indexed) in the cross-dispersion
direction, with respect to the input 2D slit image.

The coefficients for the polynomial functions are specified as a list of an
even number of lists (an even number because both the lower and upper limits of each
extraction region must be specified).  The source extraction coefficients will normally
be a list of just two lists: the coefficients for the lower limit function
and the coefficients for the upper limit function of one extraction
region.  The limits could just be constant values,
e.g. `[[324.5], [335.5]]`.  Straight but tilted lines are linear functions, e.g.
`[[324.5, 0.0137], [335.5, 0.0137]]`.

Multiple regions may be specified for either the source or background, but it is
more common to specify more than one background region.  Here
is an example for specifying two background regions:

`[[315.2, 0.0135], [320.7, 0.0135], [341.1, 0.0139], [346.8, 0.0139]]`

This is interpreted as follows:

* `[315.2, 0.0135]`: lower limit for first background region
* `[320.7, 0.0135]`: upper limit for first background region
* `[341.1, 0.0139]`: lower limit for second background region
* `[346.8, 0.0139]`: upper limit for second background region


Note that `src_coeff` and `bkg_coeff` contain floating-point
values.  For interpreting fractions of a pixel, the convention used here
is that the pixel number at the center of a pixel is a whole number.  Thus,
if a lower or upper limit is a whole number, that limit splits the pixel
in two, so the weight for that pixel will be 0.5.  To include all the
pixels between 325 and 335 inclusive, for example, the lower and upper
limits would be given as 324.5 and 335.5 respectively.

Please note that this is different from the convention used for the cross-dispersion
start/stop values, which are expected to be inclusive index values. For the example here,
for horizontal dispersion, `ystart = 325`, `ystop = 335` is equivalent
to `src_coeff = [[324.5],[335.5]]`.  To include half a pixel more at the top
and bottom of the aperture, `ystart = 324.5`, `ystop = 335.5` is equivalent
to `src_coeff = [[324],[336]]`.

The order of the polynomial is specified implicitly to be one less than the
number of coefficients. The number of coefficients for a lower or upper extraction
region limit must be at least one (i.e. zeroth-order polynomial). There is no
predefined upper limit on the number of coefficients (and hence polynomial order).
The various polynomials (lower limits, upper limits, possibly multiple regions) do
not need to have the same number of coefficients; each of the inner lists specifies
a separate polynomial. However, the independent variable (wavelength or pixel)
does need to be the same for all polynomials for a given slit.


Optimal Extraction for 2D Slit Data
-----------------------------------

Optimal extraction proceeds similarly to box extraction for 2D slit data, except that
instead of summing over an aperture defined by the reference files, a model of the point
spread function (PSF) is fit to the data at each dispersion element.  This generally provides
higher signal-to-noise for the output spectrum than box extractions and has the advantage
of ignoring missing data due to bad pixels, cosmic rays, or saturation.  Optimal extraction
also does not require a resampled spectral image as input: it can avoid the extra interpolation
by directly fitting the spatial profile along the curved trace at each dispersion element.

Optimal extraction is suited only to point sources with known source locations, for which a
high-fidelity PSF model is available.  Currently, only the MIRI LRS fixed slit exposure type
has a PSF model available in CRDS.

When optimal extraction is selected (`extraction_type = 'optimal'`), the aperture definitions in
the extraction reference file are ignored, and the following parameters
are used instead:

* `use_source_posn`: Source position is estimated from the input metadata and used to
  center the PSF model.  The recommended value is True, in order to account for spatial offsets
  within the slit; if False, or if the source position could not be estimated, the source is
  assumed to be at the center of the slit.
* `model_nod_pair`: If nod subtraction occurred prior to extraction, setting this option to
  True will allow the extraction algorithm to model a single negative trace from the nod pair
  alongside the positive trace. This can be helpful in accounting for PSF overlap between the
  positive and negative traces.  This option is ignored if no background subtraction occurred,
  or if the dither pattern was not a 2-point nod.
* `optimize_psf_location`: Since source position estimates may be slightly inaccurate,
  it may be useful to iteratively optimize the PSF location.  When this option is set to True, the
  location of the positive and negative traces (if used) are optimized with respect to the
  residuals of the scene modeled by the PSF at that location.  This option is
  strongly recommended if `model_nod_pair` is True, since the negative nod location is less
  reliably estimated than the positive trace location.
* `subtract_background`: Unlike during box extraction, the background levels can be modeled and removed
  during optimal extraction without explicitly setting a background region.  It is recommended to
  set this parameter to True if background subtraction was skipped prior to extraction. Set this
  parameter to False if a negative nod trace is present but not modeled (`model_nod_pair = False`).
* `override_psf`: If a custom flux model is required, it is possible to provide one by overriding
  the PSF model reference file. Set this parameter to the filename for a FITS file matching the
  :ref:`SpecPsfModel <psf_reffile>` format.

.. _extract-1d-for-ifu:

Extraction for 3D IFU Data
--------------------------
In IFU cube data, 1D extraction is controlled by a different set of EXTRACT1D
reference file parameters. For point source data, the extraction
aperture is centered at the RA/Dec target location indicated by the header.
If the target location is undefined in the header, then the extraction
region is the  center of the IFU cube. For extended source data, anything specified in the reference file
or step arguments will be ignored; the entire image will be extracted, and no background subtraction will be done.

For point sources, a circular extraction aperture is used, along with an optional
circular annulus for background extraction and subtraction. The size of the extraction
region and the background annulus size varies with wavelength. 
The extraction related vectors are found in the asdf extract1d reference file.
For each element in the `wavelength` vector there are three size components: `radius`, `inner_bkg`, and
`outer_bkg`. The radius vector sets the extraction size; while `inner_bkg` and `outer_bkg` specify the
limits of an annular background aperture. There are two additional vectors in the reference file, `axis_ratio`
and `axis_pa`, which are placeholders for possible future functionality.
The extraction size parameters are given in units of arcseconds and converted to units of pixels
in the extraction process. 

The region of overlap between an aperture and a pixel can be calculated by
one of three different methods, specified by the `method` parameter:  "exact"
(default), limited only by finite precision arithmetic; "center", the full value
in a pixel will be included if its center is within the aperture; or "subsample",
which means pixels will be subsampled N x N and the "center" option will be used
for each sub-pixel. When `method` is "subsample", the parameter `subpixels`
is used to set the resampling value. The default value is 10.

For IFU cubes the error information is contained entirely in the ERR array, and is not broken out into the
VAR_POISSON, VAR_RNOISE, and VAR_FLAT arrays.  As such, ``extract_1d`` only propagates this
non-differentiated error term.  Since covariance is also extremely important for undersampled IFU data
(see discussion by Law et al. 2023; AJ, 166, 45) the optional parameter `ifu_covar_scale`
will multiply all ERR arrays in the extracted spectra by a constant prefactor to account
for this covariance.  As discussed by Law et al. 2023, this prefactor provides
a reasonable first-order correction for the vast majority of use cases.  Values for the prefactor
are provided in the ``extract_1d`` parameter reference files for MIRI and NIRSpec.

.. _MIRI-MRS-1D-residual-fringe:

MIRI MRS 1D Residual Fringe Correction
--------------------------------------
For MIRI MRS IFU data there is also a correction for fringing.
As is typical for spectrometers, the MIRI MRS detectors are affected by fringes.
The primary MRS fringe, observed in all MRS bands, is caused by the etalons between the anti-reflection coating
and lower layers, encompassing the detector substrate and the infrared-active layer. Since the thickness
of the substrate is not the same in the SW and LW detectors, the fringe frequency differs in the two detectors.
Shortward of 16 microns, this fringe is produced by the anti-reflection coating and  pixel metalization etalons, whereas
longward of 16 microns it is produced by the anti-reflection coating and  bottom contact etalon, resulting in a
different fringe frequency.

The JWST pipeline contains multiple steps to mitigate the impact of fringing on science spectra and these
steps generally suffice to reduce the fringe signal to below a few percent of the target flux.

The first correction is applied by default in the :ref:`fringe <fringe_step>` step in the
:ref:`calwebb_spec2 <calwebb_spec2>` pipeline and consists of dividing the uncalibrated "rate" image
by a static fringe flat constructed from observations of a bright source that fills the entire MRS field of
view. For more details see the :ref:`fringe <fringe_step>` step.
This step generally does a good job of removing the strongest fringes from an astronomical scene, particularly
for nearly-uniform extended sources. Since the fringe signal is different for point sources, however, and varies
as a function of the location of a point source within the FOV, the static fringe flat cannot fully correct
such objects. The default high level data products will therefore still show appreciable fringes.

The pipeline also includes two optional residual fringe correction steps whose purpose is to find and remove signals
whose periodicity is consistent with known fringe frequencies (set by the optical thickness of the detectors
and dichroics) using a Lomb-Scargle periodogram. The number of fringe components to be removed is governed by
a Bayesian evidence calculation. The first of these residual fringe correction steps is a 2-D correction that
can be applied to the flux-calibrated detector data in the :ref:`residual_fringe <residual_fringe_step>` step. This step
is part of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline, but currently it is skipped by default. For more
information see :ref:`residual_fringe <residual_fringe_step>`.

The pipeline also can apply a 1-D residual fringe correction. This correction is only relevant for MIRI MRS
single band data. The parameter controlling applying the residual fringe correction is by default set to true, 
`ifu_rfcorr = True`,  in the ``extract_1d`` step.
Empirically, the 1-D correction step has been found to work better than the 2-D correction step if it is
applied to per-band spectra. If the MIRI MRS data is from multiple bands/channels the residual fringe correction
is turned off. Three additional columns are present in MIRI MRS extracted spectra: RF_FLUX, RF_SURF_BRIGHT, and 
RF_BACKGROUND. These three columns are the flux, surface brightness and background arrays with the residiual
fringe correction applied. If the data is not from a single band or the residual fringe correction fails
NaN values are reported for the arrays. 

When using the `ifu_rfcorr` option in the ``extract_1d`` step  to apply a 1-D residual fringe
correction, it is applied during the extraction of spectra from the IFU cube. The 1D residual fringe code can also
be called outside the pipeline to correct an extracted spectrum. If running outside the pipeline, the correction
works best on single-band cubes, and the channel of
the data must be given. The steps to run this correction outside the pipeline are::

  from jwst.residual_fringe.utils import fit_residual_fringes_1d as rf1d
  flux_cor = rf1d(flux, wave, channel=4)

where `flux` is the extracted spectral data, and the data are from channel 4 for this example.


Extraction for NIRISS SOSS Data
-------------------------------
For NIRISS SOSS data, spectral orders 1 and 2 overlap slightly at longer wavelengths, so a specialized extraction
algorithm known as `ATOCA (Algorithm to Treat Order ContAmination, Darveau-Bernier et al., 2022) <https://iopscience.iop.org/article/10.1088/1538-3873/ac8a77/pdf>`__ is used. This routine
constructs a linear model of each pixel on the detector and treats the underlying incident spectrum as a free variable
to simultaneously extract the cross-contaminated spectra. Using this method, the extracted spectra are accurate to
within 10ppm over the full spectral range when validated against simulations.

The algorithm uses a wavelength solution, a spectral throughput, a spectral resolution, and a spatial throughput for
both orders to determine the flux contribution from each order falling on a given pixel. Most of these references are
determined by analysis of on-sky data and supplied to the algorithm via the `pastasoss` and `spec_profile` reference
files. The exception is the `spec_kernel` reference file which supplies the convolution kernels used in the extraction,
determined from monochromatic PSFs generated by STPSF. The `pastasoss` reference file predicts the trace centroids,
taking into account the small rotations of the trace introduced by the slight visit-to-visit offsets of the GR700XD
grism in the optical path.

Having constructed a model of the intensity at each pixel using the reference file inputs, a chi-square minimization is
used to fit the pixel model to the observations, weighted by the uncertainty. However, since a solution to this system
of equations is highly degenerate, a Tikhonov regularization (Tikhonov 1963) is performed. The goal here is to find the
smoothest solution for the flux that fits the observations within the measured uncertainties.

The resulting spectral trace solutions are at a higher resolution than the observed data since an oversampled
wavelength grid is used by the ATOCA algorithm for decontamination. These results are then reconvolved onto the native
wavelength grid before the 1D spectra for each order are extracted.
