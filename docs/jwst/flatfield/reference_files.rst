Reference File 
===============
There are four reference file types for the flat_field step.  Reftype
FLAT is used for all exposure types except NIRSpec spectra.  
NIRSpec spectra use three reftypes:  FFLAT (fore optics), SFLAT (spectrograph optics), and 
DFLAT (detector).


CRDS Selection Criteria
-----------------------

For MIRI Imaging, flat-field reference files are selected based on the values of
INSTRUME, DETECTOR, FILTER, READPATT, and SUBARRAY in the science data file.

For MIRI MRS, flat-field reference files are selected based on the values of
INSTRUME, DETECTOR, BAND, READPATT, and SUBARRAY in the science data file.

For NIRCam, flat-field reference files are selected based on the values of
INSTRUME, DETECTOR, FILTER, and PUPIL in the science data file.

For NIRISS, flat-field reference files are selected based on the values of
INSTRUME, DETECTOR, and FILTER in the science data file.

For NIRSpec, flat-field reference files are selected based on the values of
INSTRUME, DETECTOR, FILTER, GRATING, and
EXP_TYPE in the science data file.


Reference File Formats for MIRI, NIRCAM, and NIRISS
---------------------------------------------------
Except for NIRSpec modes,
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
is projected).  


A single DQ_DEF extension provides the data-quality definitions for all of the 
DQ arrays, which must use the same coding scheme.  The DQ_DEF table contains 
the bit assignments used in the DQ array, and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition


Reference File Formats for NIRSpec
----------------------------------

For NIRSpec data, the flat-field reference files allow for variations in
the flat field with wavelength as well as from pixel to pixel.  There is a
separate flat-field reference file for each of three sections of the
instrument:  the fore optics (FFLAT), the spectrograph (SFLAT), and the 
detector (DFLAT).  The contents of the reference files differ from one mode 
to another (see below), but in general there may be a flat-field image and 
a 1-D array.  The image provides pixel-to-pixel values for the flat field 
that may vary slowly (or not at all) with wavelength, while the 1-D array 
is for a pixel-independent fast variation with wavelength. Details of the
file formats are given in the following sections.

If there is no significant slow variation with wavelength, the image will be a 2-D array; 
otherwise, the image will be a 3-D array, with each plane corresponding to 
a different wavelength. In the latter case, the wavelength for each plane 
will be given in a table extension called WAVELENGTH in the flat-field 
reference file.  The fast variation is given in a table extension called 
FAST_VARIATION, with column names "slit_name", "nelem", "wavelength", and 
"data" (an array of wavelength-dependent flat-field values).  Each row of 
the table contains a slit name (for fixed-slit data, otherwise "ANY"), an 
array of flat-field values, an array of the corresponding wavelengths, and 
the number of elements ("nelem") of "data" and "wavelength" that are 
populated, as the allocated array size can be larger than is needed.  
For some reference files there will not be any image array, in which case 
all the flat field information will be taken from the FAST_VARIATION table.  

The SCI extension of the reference file may contain NaNs.  If so, the
flat_field step will replace these values with 1 and will flag the
corresponding pixel in the DQ extension with NO_FLAT_FIELD.  The WAVELENGTH
extension is not expected to contain NaNs.

For the detector section, there is only one flat-field reference file for
each detector.  For the fore optics and the spectrograph sections, however,
there are different flat fields for fixed-slit data, IFU data, and for
multi-object spectroscopic data.  Here is a summary of the contents of these
files.

For the fore optics, the flat field for fixed-slit data contains just a
FAST_VARIATION table (i.e. there is no image).  This table has five rows,
one for each of the fixed slits.  The flat field for IFU data also contains
just a FAST_VARIATION table, but it has only one row with the value "ANY"
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
(i.e. the flat field at multiple wavelengths), a corresponding
WAVELENGTH table, and a FAST_VARIATION table with one row.

As just described, there are 3 types of reference files for NIRSpec (FFLAT, 
SFLAT, and DFLAT), and within each of these types, there are several formats, 
which are now described.


Fore Optics (FFLAT)
:::::::::::::::::::
There are 3 types of FFLAT reference files: fixed slit, msa spec, and IFU. For each type
the primary data array is assumed to be empty.


*Fixed Slit*
~~~~~~~~~~~~
The fixed slit references files have EXP_TYPE=NRS_FIXEDSLIT, and have a single BINTABLE
extension, labeled FAST_VARIATION. 

The table contains four columns:

* slit_name: string, name of slit
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength

The number of rows in the table is given by NAXIS2, and each row corresponds to a separate slit.


*MSA Spec*
~~~~~~~~~~
The MSA Spec references files have EXP_TYPE=NRS_MSASPEC, and contain data pertaining
to each of the 4 quadrants.  For each quadrant, there are 3 IMAGE extensions, a BINTABLE extension 
labeled WAVELENGTH, and a BINTABLE extension labeled FAST_VARIATION.  The file also contains 
one BINTABLE labeled DQ_DEF.

The IMAGE extensions have the following characteristics:

=======   =====  =====================  =========
EXTNAME   NAXIS  Dimensions             Data type
=======   =====  =====================  =========
SCI       3      ncols x nrows x nelem  float
ERR       3      ncols x nrows x nelem  float
DQ        3      ncols x nrows x nelem  integer
=======   =====  =====================  =========

For all 3 of these extensions, the EXTVER keyword indicates the quadrant number, 1 to 4.
Each plane of the SCI array gives the flat_field value for every pixel in the quadrant for
the corresponding wavelength, which is specified in the WAVELENGTH table.



The WAVELENGTH table contains a single column:

* wavelength: float 1-D array, values of wavelength

Each of these wavelength values corresponds to a single plane of the IMAGE arrays.


The FAST_VARIATION table contains four columns:

* slit_name: the string "ANY"
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength


The flat field values in this table are used to account for a wavelength-dependence on a much
finer scale than given by the values in the SCI array.  There is a single row in this table, 
as the same wavelength-dependent value is applied to all pixels in the quadrant.

 
The DQ_DEF table contains the bit assignments used in the DQ array, and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition


*IFU*
~~~~~
The IFU reference files have EXP_TYPE=NRS_IFU.  These have one extensions,
a BINTABLE extension labeled FAST_VARIATION. 

The FAST_VARIATION table contains four columns:

* slit_name: the string "ANY"
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength

The flat field values in this table are used to account for a 
wavelength-dependence on a much finer scale than given by the values in the SCI 
array. For each pixel in the science data, the wavelength of the light that fell
on that pixel will be determined by using the WCS interface. The flat-field 
value for that pixel will then be obtained by interpolating within the 
wavelength and data arrays from the FAST_VARIATION table.

The DQ_DEF table contains the bit assignments used in the DQ arrays. The table contains the 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition


Spectrograph (SFLAT)
::::::::::::::::::::

There are 3 types of SFLAT reference files: fixed slit, msa spec, and IFU. For each type
the primary data array is assumed to be empty.


*Fixed Slit*
~~~~~~~~~~~~
The fixed slit references files have EXP_TYPE=NRS_FIXEDSLIT, and have a BINTABLE
extension labeled FAST_VARIATION. The table contains four columns:

* slit_name: string, name of slit
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength

The number of rows in the table is given by NAXIS2, and each row corresponds to a separate slit.


*MSA Spec*
~~~~~~~~~~
The MSA Spec references files have EXP_TYPE=NRS_MSASPEC. There are 3 IMAGE extensions, a BINTABLE extension 
labeled WAVELENGTH, a BINTABLE extension labeled FAST_VARIATION, and a BINTABLE labeled DQ_DEF.

The IMAGE extensions have the following characteristics:

=======   =====  ====================  =========
EXTNAME   NAXIS  Dimensions             Data type
=======   =====  ====================  =========
SCI       3      ncols x nrows x n_wl  float
ERR       3      ncols x nrows x n_wl  float
DQ        3      ncols x nrows x n_wl  integer
=======   =====  ====================  =========

The keyword NAXIS3 in these extensions specifies the number n_wl of monochromatic slices, each of which
gives the flat_field value for every pixel for the corresponding wavelength, which is 
specified in the WAVELENGTH table.


The WAVELENGTH table contains a single column:

* wavelength: float 1-D array, values of wavelength

Each of these wavelength values corresponds to a single plane of the IMAGE arrays.


The FAST_VARIATION table contains four columns:

* slit_name: the string "ANY"
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength

The flat field values in this table are used to account for a wavelength-dependence on a much
finer scale than given by the values in the SCI array.  For each pixel in the science data, 
the wavelength of the light that fell on that pixel will be determined by using the WCS
interface.  The flat-field value for that pixel will then be obtained by
interpolating within the wavelength and data arrays from the FAST_VARIATION
table.

 
The DQ_DEF table contains the bit assignments used in the DQ array, and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition


Detector (DFLAT)
::::::::::::::::

There is only one type of DFLAT reference file, and it contains 3 IMAGE extensions, a BINTABLE extension 
labeled WAVELENGTH, a BINTABLE extension labeled FAST_VARIATION, and a BINTABLE labeled DQ_DEF.

The IMAGE extensions have the following characteristics:

=======   =====  ====================  =========
EXTNAME   NAXIS  Dimensions            Data type
=======   =====  ====================  =========
SCI       3      ncols x nrows x n_wl  float
ERR       3      ncols x nrows         float
DQ        3      ncols x nrows         integer
=======   =====  ====================  =========


The keyword NAXIS3 in the SCI IMAGE extension specifies the number n_wl of monochromatic slices, 
each of which gives the flat_field value for every pixel for the corresponding wavelength, which is 
specified in the WAVELENGTH table.

The WAVELENGTH table contains a single column:

* wavelength: float 1-D array, values of wavelength

Each of these wavelength values corresponds to a single plane of the SCI IMAGE array.

The FAST_VARIATION table contains four columns:

* slit_name: the string "ANY"
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength


The flat field values in this table are used to account for a wavelength-dependence on a much
finer scale than given by the values in the SCI array.  There is a single row in this table, 
as the same wavelength-dependent value is applied to all pixels.

The DQ_DEF table contains the bit assignments used in the DQ array, and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition
