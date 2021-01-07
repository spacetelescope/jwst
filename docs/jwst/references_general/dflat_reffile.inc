.. _dflat_reffile:

DFLAT Reference File
--------------------

:REFTYPE: DFLAT
:Data model: `~jwst.datamodels.NirspecFlatModel`

.. include:: ../references_general/dflat_selection.inc

There is one type of DFLAT reference file, containing 6 extensions with the following
characteristics:

============== ======== ===== ===================== =========
EXTNAME        XTENSION NAXIS Dimensions            Data type
============== ======== ===== ===================== =========
SCI            IMAGE      3   ncols x nrows x n_wl  float
ERR            IMAGE      2   ncols x nrows         float
DQ             IMAGE      2   ncols x nrows         integer
WAVELENGTH     BINTABLE   2   TFIELDS = 1           N/A
FAST_VARIATION BINTABLE   2   TFIELDS = 4           N/A
DQ_DEF         BINTABLE   2   TFIELDS = 4           N/A
============== ======== ===== ===================== =========

.. include:: ../includes/dq_def.inc

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
this table, which contains 1-D arrays of wavelength and flat-field values.
The same wavelength-dependent value is applied to all pixels in a quadrant.

