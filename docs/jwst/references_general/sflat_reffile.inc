.. _sflat_reffile:

SFLAT Reference File
--------------------

:REFTYPE: SFLAT
:Data model: `~jwst.datamodels.NirspecFlatModel`

There are 3 types of NIRSpec SFLAT reference files: fixed slit, MSA spec, and IFU.
For each type the primary HDU does not contain a data array.

.. include:: ../references_general/sflat_selection.inc

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

============== ======== ===== ===================== =========
EXTNAME        XTENSION NAXIS Dimensions            Data type
============== ======== ===== ===================== =========
SCI            IMAGE      3   ncols x nrows x n_wl  float
ERR            IMAGE      3   ncols x nrows x n_wl  float
DQ             IMAGE      3   ncols x nrows x n_wl  integer
WAVELENGTH     BINTABLE   2   TFIELDS = 1           N/A
FAST_VARIATION BINTABLE   2   TFIELDS = 4           N/A
DQ_DEF         BINTABLE   2   TFIELDS = 4           N/A
============== ======== ===== ===================== =========

.. include:: ../includes/dq_def.inc

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
science data, the wavelength of the light that fell on that pixel will be read from the
WAVELENGTH array in the science exposure (if that array is absent, it will be computed
using the WCS interface).
The flat-field value for that pixel will then be obtained by
interpolating within the wavelength and data arrays from the FAST_VARIATION
table.
