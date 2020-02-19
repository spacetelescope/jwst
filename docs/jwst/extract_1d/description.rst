Description
===========
The extract_1d step extracts a 1-d signal from a 2-d or 3-d dataset and
writes a spectrum to a product.  This works on fixed-slit data (NIRSpec
data through any one or more of the fixed slits, MIRI LRS data through
the slit or in the slitless region, and NIRISS slitless data) as well as
IFU data and NIRSpec MOS (micro-shutter array) data.

For GRISM data (NIS_WFSS or NRC_WFSS), no reference file is used.
The extraction region is taken to be the full size of the input subarray
or cutout, or it could be restricted to the region within which the
world coordinate system is defined.  The dispersion direction is the one
along which the wavelengths change more rapidly.

For IFU data, the extraction options differ depending on
whether the target is a point source or an extended source.  For a point
source, the spectrum will be extracted using circular aperture photometry,
optionally including background subtraction using a circular annulus.
For an extended source, rectangular aperture photometry will be used, with
the entire image being extracted, and no background subtraction, regardless
of what was specified in the reference file or command-line arguments.
For either point source or extended, the photometry makes use of
astropy photutils.
The region of overlap between an aperture and a pixel can be calculated by
one of three different methods:  "exact", limited only by finite precision
arithmetic; "center", i.e. the full value in a pixel will be included if its
center is within the aperture; or "subsample", which means pixels will be
subsampled N x N, and the "center" option will be used for each sub-pixel.


Input
=====
Level 2-b countrate data, or level-3 data.  The format should be a
CubeModel, a SlitModel, an IFUCubeModel, an ImageModel, a DrizProductModel,
a MultiSlitModel, a MultiProductModel, or a ModelContainer.
The SCI extensions should
have keyword SLTNAME to specify which slit was extracted, though if there
is only one slit (e.g. full-frame data), the slit name can be taken from
the JSON reference file instead.

For almost all purposes, the photom step should have been run before running
extract_1d.  (If photom has not been run, a warning will be logged and the
output of extract_1d will be in units of count rate.)  The photom step
reads data in units of count rate and converts to units of either surface
brightness (megajanskys per steradian) or, for point sources observed with
NIRSpec or NIRISS SOSS only, units of flux density (megajanskys).

Output
======
The output will be in MultiSpecModel format; for each input slit there will
be an output table extension with the name EXTRACT1D.  This extension will
have columns WAVELENGTH, FLUX, ERROR, SURF_BRIGHT, SB_ERROR, DQ,
BACKGROUND, BERROR and NPIXELS.

Some metadata will be written to the table header, mostly copied from the
input header.

WAVELENGTH was copied from the wavelength attribute of the input 2-D data,
if that attribute exists and was populated, or it was calculated from the
WCS otherwise.
FLUX is the flux density in janskys; see keyword TUNIT2 if the data are
in a FITS BINTABLE.  ERROR is the error estimate for FLUX, and it has the
same units as FLUX.
SURF_BRIGHT is the surface brightness in MJy / sr, except that for point
sources observed with NIRSpec or NIRISS SOSS, SURF_BRIGHT will be set to
zero.  SB_ERROR is the error estimate for SURF_BRIGHT.
While it was expected that a user would make use of the FLUX column for
point-source data or the SURF_BRIGHT column for an extended source,
both columns will be populated, except for NIRSpec or NIRISS SOSS for
point sources.
The extract_1d step collapses the input data from 2-D to 1-D by summing
one or more rows (or columns, depending on the dispersion direction).
A background may optionally be subtracted by the extract_1d step, but
there are also other options for background subtraction prior to extract_1d.
For the case that the input data are in units of MJy / sr, the SURF_BRIGHT
and BACKGROUND columns will be
populated by dividing the sum by the number of pixels (see the NPIXELS column,
described below) that were added together; the FLUX column will be populated
by multiplying the sum by the solid angle of a pixel, and also multiplying
by 10^6 to convert from MJy to Jy.
For the case that the input data are in units of MJy (i.e. point sources,
NIRSpec or NIRISS SOSS), the SURF_BRIGHT column will be set to zero; the
FLUX column will just be multiplied by 10^6; the BACKGROUND column will be
divided by NPIXELS and by the solid angle of a pixel to convert to surface
brightness (MJy / sr).

NPIXELS is the number of pixels that were added together for the source
extraction region.  Note that this is not necessarily a constant, and
the value is not necessarily an integer (the data type is float).
BACKGROUND is the measured background, scaled to the extraction width used
for FLUX and SURF_BRIGHT.  BACKGROUND will be zero if no background was
subtracted in the extract_1d step.
ERROR, SB_ERROR, BERROR, and DQ are not populated with useful values yet.
