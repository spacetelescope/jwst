Reference File
==============
The flat-field correction step uses a FLAT reference file.

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
Except for NIRSpec MOS mode (see below),
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

Reference File Format for NIRSpec MOS data
------------------------------------------

For the NIRSpec multi-object spectroscopy (MOS) mode, the flat field
reference file contains three IMAGE extensions, SCI, WAVELENGTH, and DQ,
and one BINTABLE extension, DQ_DEF.  The detector response varies with
wavelength, and light of different wavelengths can fall on each pixel, so
for a given science exposure, the flat field that should be applied to any
given pixel will be obtained by interpolating within a set of flat-field
reference images that cover different ranges of wavelength.  The SCI,
WAVELENGTH, and DQ extensions in the reference file are 3-D images, all
the same shape.  The SCI array is a stack of several (perhaps a few dozen)
full-frame flat-field images.  The WAVELENGTH array gives the wavelength
of the light that illuminated the detector for each pixel of each flat-field
image.  The DQ array flags pixels for which the data value is not valid
(NO_FLAT_FIELD).

The SCI extension of the reference file may contain NaNs.  If so, the
flat_field step will replace these values with 1 and will flag the
corresponding pixel in the DQ extension with NO_FLAT_FIELD.  The WAVELENGTH
extension is not expected to contain NaNs.
