.. _fflat_reffile:

FFLAT Reference File
--------------------

:REFTYPE: FFLAT

There are 3 forms of NIRSpec FFLAT reference files: fixed slit, MSA spec, and IFU.
For each type the primary HDU does not contain a data array.

.. include:: ../references_general/fflat_selection.inc

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

============== ======== ===== ===================== =========
EXTNAME        XTENSION NAXIS Dimensions            Data type
============== ======== ===== ===================== =========
SCI            IMAGE      3   ncols x nrows x nelem float
ERR            IMAGE      3   ncols x nrows x nelem float
DQ             IMAGE      3   ncols x nrows x nelem integer
WAVELENGTH     BINTABLE   2   TFIELDS = 1           N/A
FAST_VARIATION BINTABLE   2   TFIELDS = 4           N/A
DQ_DEF         BINTABLE   2   TFIELDS = 4           N/A
============== ======== ===== ===================== =========

.. include:: ../includes/dq_def.inc

For the 5 extensions that appear multiple times, the EXTVER keyword indicates the
quadrant number, 1 to 4.
Each plane of the SCI array gives the throughput value for every shutter in the
MSA quadrant for the corresponding wavelength, which is specified in the
WAVELENGTH table. These wavelength-dependent values are combined with the
FAST_VARIATION array, and are then applied to the science spectrum based on the
wavelength of each pixel.

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
this table, which contains 1-D arrays of wavelength and flat-field values.
The same wavelength-dependent value is applied to all pixels in a quadrant.

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

For each pixel in the science data, the wavelength of the light that fell
on that pixel will be determined from the WAVELENGTH array in the science exposure
(in the absence of that array, it will be computed using the WCS interface).
The flat-field value for that pixel will then be obtained by interpolating within the 
wavelength and data arrays from the FAST_VARIATION table.
