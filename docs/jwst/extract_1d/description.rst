Description
===========
The ``extract_1d`` step extracts a 1D signal from a 2D or 3D dataset and
writes a spectrum to a product.  This works on all JWST spectroscopic
modes, including MIRI LRS (slit and slitless) and MRS, NIRCam WFSS and
TSGRISM, NIRISS WFSS and SOSS, and NIRSpec fixed-slit, IFU, and MOS.

An EXTRACT1D reference file is used for most modes to specify the location and
size of the target and background extraction apertures.
The EXTRACT1D reference file is not used for Wide-Field Slitless Spectroscopy data
(NIS_WFSS or NRC_WFSS). The extraction region is instead taken to be the full size
of the input subarray or cutout, or restricted to the region within which the world
coordinate system (WCS) is defined.

For IFU data, the extraction options differ depending on
whether the target is a point source or an extended source.  For a point
source, the spectrum is extracted using circular aperture photometry,
optionally including background subtraction using a circular annulus.
For an extended source, rectangular aperture photometry is used, with
the entire image being extracted, and no background subtraction, regardless
of what was specified in the reference file or command-line arguments.
For either point source or extended, the photometry makes use of astropy photutils.
The region of overlap between an aperture and a pixel can be calculated by
one of three different methods:  "exact", limited only by finite precision
arithmetic; "center", i.e. the full value in a pixel will be included if its
center is within the aperture; or "subsample", which means pixels will be
subsampled N x N, and the "center" option will be used for each sub-pixel.

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
The output will be in MultiSpecModel format. For each input slit there will
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
