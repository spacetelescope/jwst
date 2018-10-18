Reference Files
===============
The ``flat_field`` step uses four different types of reference files, depending on the
type of data being processed. Most cases just use the FLAT reference file, while NIRSpec
spectroscopic exposures use the three reference files FFLAT (fore optics),
SFLAT (spectrograph optics), and DFLAT (detector).

FLAT Reference File
-------------------

:REFTYPE: FLAT
:Data model: `~jwst.datamodels.FlatModel`

The FLAT reference file contains pixel-by-pixel detector response values.
It is used for all instrument modes except the NIRSpec spectroscopic modes.

.. include:: flat_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for FLAT
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in FLAT reference files,
because they are used as CRDS selectors
(see :ref:`flat_selectors`):

=========  ==============================  =============================
Keyword    Data Model Name                 Instruments
=========  ==============================  =============================
DETECTOR   model.meta.instrument.detector  All
EXP_TYPE   model.meta.exposure.type        FGS, NIRSpec
FILTER     model.meta.instrument.filter    MIRI, NIRCam, NIRISS, NIRSpec
PUPIL      model.meta.instrument.pupil     NIRCam, NIRISS
BAND       model.meta.instrument.band      MIRI
READPATT   model.meta.exposure.readpatt    MIRI
SUBARRAY   model.meta.subarray.name        MIRI
GRATING    model.meta.instrument.grating   NIRSpec
=========  ==============================  =============================

Reference File Format
+++++++++++++++++++++
FLAT reference files are FITS format, with 3 IMAGE extensions
and 1 BINTABLE extension. The FITS primary data array is assumed to be empty.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
ERR      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   integer
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

The ``DQ_DEF`` table extension lists the bit assignments for the flag conditions
used in the DQ array.

.. include:: ../includes/dq_def.rst

For application to imaging data, the FITS file contains a single set of SCI,
ERR, DQ, and DQ_DEF extensions.  Image dimensions should be 2048x2048 for the
NIR detectors and 1032x1024 for MIRI (i.e. they include reference pixels),
unless data were taken in subarray mode.

Reference Files for NIRSpec Spectroscopy
===================================================
For NIRSpec spectroscopic data, the flat-field reference files allow for variations in
the flat field with wavelength, as well as from pixel to pixel.  There is a
separate flat-field reference file for each of three sections of the
instrument:  the fore optics (FFLAT), the spectrograph (SFLAT), and the
detector (DFLAT).  The contents of the reference files differ from one mode
to another (see below), but in general they may contain a flat-field image and
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
populated, because the allocated array size can be larger than needed.
For some reference files there will not be any image array, in which case
all the flat field information will be taken from the FAST_VARIATION table.

The SCI extension of the reference files may contain NaNs.  If so, the
``flat_field`` step will replace these values with 1 and will flag the
corresponding pixel in the DQ extension with NO_FLAT_FIELD. The WAVELENGTH
extension is not expected to contain NaNs.

For the detector section, there is only one flat-field reference file for
each detector.  For the fore optics and the spectrograph sections, however,
there are different flat fields for fixed-slit data, IFU data, and for
multi-object spectroscopic data.  Here is a summary of the contents of these
files.

For the fore optics (FFLAT), the flat field for fixed-slit data contains just a
FAST_VARIATION table (i.e. there is no image).  This table has five rows,
one for each of the fixed slits.  The FFLAT for IFU data also contains
just a FAST_VARIATION table, but it has only one row with the value "ANY"
in the "slit_name" column.  For multi-object spectroscopic data, the FFLAT
contains four sets of images (one for each MSA quadrant), WAVELENGTH
tables, and FAST_VARIATION tables.  The images are unique to the FFLATs,
however.  The image "pixels" correspond to micro-shutter array slits, rather
than to detector pixels.  The array size is 365 rows
by 171 columns, and there are multiple planes to handle the slow variation
of flat field with wavelength.

For the spectrograph optics (SFLAT), the flat-field files have nearly the same
format for fixed-slit data, IFU, and multi-object data.  The difference is
that for fixed-slit and IFU data, the image is just a single plane,
i.e. the only variation with wavelength is in the FAST_VARIATION table,
while there are multiple planes in the image for multi-object spectroscopic
data (and therefore there is also a corresponding WAVELENGTH table, with
one row for each plane of the image).

For the detector section, the DFLAT file contains a 3-D image
(i.e. the flat field at multiple wavelengths), a corresponding
WAVELENGTH table, and a FAST_VARIATION table with one row.

As just described, there are 3 types of reference files for NIRSpec (FFLAT, 
SFLAT, and DFLAT), and within each of these types, there are several formats, 
which are now described.

FFLAT Reference File
--------------------

:REFTYPE: FFLAT

There are 3 forms of FFLAT reference files: fixed slit, MSA spec, and IFU. For each type
the primary data array is assumed to be empty.

.. include:: fflat_selection.rst

*Fixed Slit*
++++++++++++
:Data model: `~jwst.datamodels.NirspecFlatModel`

The fixed slit FFLAT files have EXP_TYPE=NRS_FIXEDSLIT, and have a single BINTABLE
extension, labeled FAST_VARIATION. 

The table contains four columns:

* slit_name: string, name of slit
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength

The number of rows in the table is given by NAXIS2, and each row corresponds to a
separate slit.

*MSA Spec*
++++++++++++
:Data model: `~jwst.datamodels.NirspecQuadFlatModel`

The MSA Spec FFLAT files have EXP_TYPE=NRS_MSASPEC, and contain data pertaining
to each of the 4 quadrants.  For each quadrant, there is a set of 5 extensions -
SCI, ERR, DQ, WAVELENGTH, and FAST_VARIATION.
The file also contains a single DQ_DEF extension.

The extensions have the following characteristics:

============   ======== ===== ===================== =========
EXTNAME        XTENSION NAXIS Dimensions            Data type
=============  ======== ===== ===================== =========
SCI            IMAGE      3   ncols x nrows x nelem float
ERR            IMAGE      3   ncols x nrows x nelem float
DQ             IMAGE      3   ncols x nrows x nelem integer
WAVELENGTH     BINTABLE   2   TFIELDS = 1           N/A
FAST_VARIATION BINTABLE   2   TFIELDS = 4           N/A
DQ_DEF         BINTABLE   2   TFIELDS = 4           N/A
============== ======== ===== ===================== =========

.. include:: ../includes/dq_def.rst

For the 5 extensions that appear multiple times, the EXTVER keyword indicates the
quadrant number, 1 to 4.
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

The flat field values in this table are used to account for a wavelength-dependence on a
much finer scale than given by the values in the SCI array.  There is a single row in
this table, as the same wavelength-dependent value is applied to all pixels in the
quadrant.

*IFU*
+++++
:Data model: `~jwst.datamodels.NirspecFlatModel`

The IFU FFLAT files have EXP_TYPE=NRS_IFU.  These have one extension,
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

SFLAT Reference File
--------------------

:REFTYPE: SFLAT
:Data model: `~jwst.datamodels.NirspecFlatModel`

There are 3 types of SFLAT reference files: fixed slit, MSA spec, and IFU. For each type
the primary data array is assumed to be empty.

.. include:: sflat_selection.rst

*Fixed Slit*
++++++++++++
The fixed slit references files have EXP_TYPE=NRS_FIXEDSLIT, and have a BINTABLE
extension labeled FAST_VARIATION. The table contains four columns:

* slit_name: string, name of slit
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength

The number of rows in the table is given by NAXIS2, and each row corresponds to a
separate slit.

*MSA Spec*
++++++++++
The MSA Spec SFLAT files have EXP_TYPE=NRS_MSASPEC.
They contain 6 extensions, with the following characteristics:

============   ======== ===== ===================== =========
EXTNAME        XTENSION NAXIS Dimensions            Data type
=============  ======== ===== ===================== =========
SCI            IMAGE      3   ncols x nrows x n_wl  float
ERR            IMAGE      3   ncols x nrows x n_wl  float
DQ             IMAGE      3   ncols x nrows x n_wl  integer
WAVELENGTH     BINTABLE   2   TFIELDS = 1           N/A
FAST_VARIATION BINTABLE   2   TFIELDS = 4           N/A
DQ_DEF         BINTABLE   2   TFIELDS = 4           N/A
============== ======== ===== ===================== =========

.. include:: ../includes/dq_def.rst

The keyword NAXIS3 in the 3 IMAGE extensions specifies the number, n_wl, of monochromatic
slices, each of which gives the flat_field value for every pixel for the corresponding
wavelength, which is specified in the WAVELENGTH table.

The WAVELENGTH table contains a single column:

* wavelength: float 1-D array, values of wavelength

Each of these wavelength values corresponds to a single plane of the IMAGE arrays.

The FAST_VARIATION table contains four columns:

* slit_name: the string "ANY"
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength

The flat field values in this table are used to account for a wavelength-dependence on a
much finer scale than given by the values in the SCI array.  For each pixel in the
science data, the wavelength of the light that fell on that pixel will be determined by
using the WCS interface.  The flat-field value for that pixel will then be obtained by
interpolating within the wavelength and data arrays from the FAST_VARIATION
table.

DFLAT Reference File
--------------------

:REFTYPE: DFLAT
:Data model: `~jwst.datamodels.NirspecFlatModel`

.. include:: dflat_selection.rst

There is only one type of DFLAT reference file, containing 6 extensions with the following
characteristics:

============   ======== ===== ===================== =========
EXTNAME        XTENSION NAXIS Dimensions            Data type
=============  ======== ===== ===================== =========
SCI            IMAGE      3   ncols x nrows x n_wl  float
ERR            IMAGE      2   ncols x nrows         float
DQ             IMAGE      2   ncols x nrows         integer
WAVELENGTH     BINTABLE   2   TFIELDS = 1           N/A
FAST_VARIATION BINTABLE   2   TFIELDS = 4           N/A
DQ_DEF         BINTABLE   2   TFIELDS = 4           N/A
============== ======== ===== ===================== =========

.. include:: ../includes/dq_def.rst

The keyword NAXIS3 in the SCI extension specifies the number, n_wl, of monochromatic
slices, each of which gives the flat_field value for every pixel for the corresponding
wavelength, which is specified in the WAVELENGTH table.

The WAVELENGTH table contains a single column:

* wavelength: float 1-D array, values of wavelength

Each of these wavelength values corresponds to a single plane of the SCI IMAGE array.

The FAST_VARIATION table contains four columns:

* slit_name: the string "ANY"
* nelem: integer, maximum number of wavelengths
* wavelength: float 1-D array, values of wavelength
* data: float 1-D array, flat field values for each wavelength

The flat field values in this table are used to account for a wavelength-dependence on a
much finer scale than given by the values in the SCI array.  There is a single row in
this table, because the same wavelength-dependent value is applied to all pixels.

