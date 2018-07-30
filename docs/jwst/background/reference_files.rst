Reference Files
===============
The background image subtraction step uses reference files only when
processing Wide-Field Slitless Spectroscopy (WFSS) exposures. Two reference
files are used for WFSS mode.

WFSS Background reference file
------------------------------

:REFTYPE: WFSSBKG
:Data model: `WfssBkgModel`

The WFSS background reference file contains a "master" image of the
dispersed background produced by a particular filter+grism combination.

CRDS Selection Criteria
+++++++++++++++++++++++
WFSSBKG reference files are selected by:

  INSTRUME, DETECTOR, EXP_TYPE, FILTER, and PUPIL

Required Keywords
+++++++++++++++++
The following table lists the keywords that are required to be present in
a WFSSBKG reference file. An asterisk following a keyword name indicates a
standard keyword that is required in all reference files, regardless of
type.

=========  ========================
Keyword    Model Name
=========  ========================
AUTHOR*    meta.author
DATAMODL*  meta.model_type
DATE*      meta.date
DESCRIP*   meta.description
DETECTOR   meta.instrument.detector
EXP_TYPE   meta.exposure.type
FILENAME*  meta.filename
FILTER     meta.instrument.filter
INSTRUME*  meta.instrument.name
PEDIGREE*  meta.pedigree
PUPIL      meta.instrument.pupil
REFTYPE*   meta.reftype
TELESCOP*  meta.telescope
USEAFTER*  meta.useafter
=========  ========================

Reference File Format
+++++++++++++++++++++
WFSSBKG reference files are FITS files with 3 IMAGE extensions and
1 BINTABLE extension. The FITS primary data array is assumed to be empty.
The characteristics of the FITS extensions are as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
ERR      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   integer
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

The DQ_DEF extension contains the bit assignments used in the DQ array.
It contains the following 4 columns:

===========  =======  ===============================================
TTYPE        TFORM    Description
===========  =======  ===============================================
BIT          integer  The bit number, starting at zero
VALUE        integer  The equivalent base-10 value of BIT
NAME         string   The mnemonic name of the data quality condition
DESCRIPTION  string   A description of the data quality condition
===========  =======  ===============================================

Wavelength Range reference file
-------------------------------

:REFTYPE: WAVELENGTHRANGE
:Data model: `WavelengthrangeModel`

The wavelength range reference file contains information about the range of
wavelengths in the exposure. It is used, together with a source catalog,
to create a mask giving the locations of source spectra in the target image
and hence where the background regions are.

CRDS Selection Criteria
+++++++++++++++++++++++
Wavelengthrange reference files are selected by:

  INSTRUME, EXP_TYPE, PUPIL (NIRCam only), and MODULE (NIRCam only)

