Reference File
==============
The flat-field correction step uses a FLAT reference file.
For NIRSpec data, however, there are three flat-field reference files,
FFLAT (fore optics), SFLAT (spectrograph optics), and
DFLAT (detector).

CRDS Selection Criteria
-----------------------
Flat-field reference files are selected by the following criteria:

- MIRI Imager: Match INSTRUME, DETECTOR, FILTER, READPATT, and
  SUBARRAY of the science data file.

- MIRI MRS: Match INSTRUME, DETECTOR, BAND, READPATT, and
  SUBARRAY of the science data file.

- NIRCam: Match INSTRUME, DETECTOR, FILTER, and PUPIL.

- NIRISS: Match INSTRUME, DETECTOR, and FILTER.

- NIRSpec: Match INSTRUME, DETECTOR, FILTER, GRATING, and
  EXP_TYPE.

Reference File Format
---------------------
Except for NIRSpec data (see below),
flat-field reference files are FITS format with 3 IMAGE extensions and 1
BINTABLE extension. The primary data array is assumed to be empty. The 3
IMAGE extensions have the following characteristics:

========  =====  =============  =========
EXTNAME   NAXIS  Dimensions     Data type
========  =====  =============  =========
SCI       2      ncols x nrows  float
ERR       2      ncols x nrows  float
DQ        2      ncols x nrows  integer
========  =====  =============  =========

The BINTABLE extension uses ``EXTNAME=DQ_DEF`` and contains the bit assignments
of the conditions flagged in the DQ array.

For application to imaging data, the FITS file contains a single set of SCI,
ERR, DQ, and DQ_DEF extensions.  Image dimensions should be 2048x2048 for the
NIR detectors and 1032 x 1024 for MIRI, unless data were taken in subarray
mode.

For slit spectroscopy, a set of SCI, ERR and DQ extensions can be provided
for each aperture (identified by the detector subarray onto which the spectrum
is projected).  A single DQ_DEF extension provides the data-quality definitions
for all of the DQ arrays, which must use the same coding scheme.

Reference File Format for NIRSpec data
------------------------------------------
For NIRSpec data, the flat-field reference files allow for variations in
the flat field with wavelength as well as from pixel to pixel.  There is a
separate flat-field reference file for each of three sections of the
instrument:  the fore optics, the spectrograph, and the detector.  The
contents of the reference files differ from one mode to another (see below),
but in general there may be a flat-field image and a 1-D array.  The image
provides pixel-to-pixel values for the flat field that may vary slowly (or
not at all) with wavelength, while the 1-D array is for a pixel-independent
fast variation with wavelength.  If there is no significant slow variation
with wavelength, the image will be a 2-D array; otherwise, the image will
be a 3-D array, with each plane corresponding to a different wavelength.
In the latter case, the wavelength for each plane will be given in a table
extension called WAVELENGTH in the flat-field reference file.  The fast
variation is given in a table extension called FAST_VARIATION, with column
names "slit_name", "nelem", "wavelength", and "data" (an array of
wavelength-dependent flat-field values).  Each row of the table contains
a slit name (for fixed-slit data, otherwise "ANY"), an array of flat-field
values, an array of the corresponding wavelengths, and the number of
elements ("nelem") of "data" and "wavelength" that are populated, as the
allocated array size can be larger than is needed.  For some reference
files there will not be any image array, in which case all the flat field
information will be taken from the FAST_VARIATION table.

For the detector section, there is only one flat-field reference file for
each detector.  For the fore optics and the spectrograph sections, however,
there are different flat fields for fixed-slit data, IFU data, and for
multi-object spectroscopic data.  Here is a summary of the contents of these
files.

For the fore optics, the flat field for fixed-slit data contains just a
FAST_VARIATION table (i.e. there is no image).  This table has five rows,
one for each of the fixed slits.  The flat field for IFU data also contains
just a FAST_VARIATION table, but it has only one row (with the value "ANY"
in the "slit_name" column.  For multi-object spectroscopic data, the flat
field contains four sets (one for each MSA quadrant) of images, WAVELENGTH
tables, and FAST_VARIATION tables.  The images are unique to the fore
optics flat fields, however.  The image "pixels" correspond to micro-shutter
array slits, rather than to detector pixels.  The array size is 365 rows
by 171 columns, and there are multiple planes to handle the slow variation
of flat field with wavelength.

For the spectrograph optics, the flat-field files have nearly the same
format for fixed-slit data, IFU, and multi-object data.  The difference is
that for fixed-slit and IFU data, the image is just a single plane,
i.e. the only variation with wavelength is in the FAST_VARIATION table,
while there are multiple planes in the image for multi-object spectroscopic
data (and therefore there is also a corresponding WAVELENGTH table, with
one row for each plane of the image).

For the detector section, the flat field file contains a 3-D image
(i.e. the the flat field at multiple wavelengths), a corresponding
WAVELENGTH table, and a FAST_VARIATION table with one row.
