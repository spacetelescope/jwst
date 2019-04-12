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
no background subtraction.  The photometry makes use of astropy photutils.
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

Output
======
The output will be in MultiSpecModel format; for each input slit there will
be an output table extension with the name EXTRACT1D.  This extension will
have columns WAVELENGTH, FLUX, ERROR, DQ, NET, NERROR, BACKGROUND, BERROR
and NPIXELS.
WAVELENGTH was copied from the wavelength attribute of the input 2-D data,
if that attribute exists and was populated, or it was calculated from the
WCS otherwise.
NET is the count rate minus background, in counts/second (per pixel in the
dispersion direction), obtained by summing along the direction perpendicular
to the dispersion.  Currently only a simple summation is done, with
no weighting.
NPIXELS is the number of pixels that were added together for the source
extraction region.
BACKGROUND is the measured background, scaled to the extraction width used
for the NET.  BACKGROUND will be zero if no background was subtracted.
FLUX will be computed from NET if there is a RELSENS table
for the input slit; otherwise, FLUX will be zero.
ERROR, DQ, NERROR, and BERROR are not populated with useful values yet.
