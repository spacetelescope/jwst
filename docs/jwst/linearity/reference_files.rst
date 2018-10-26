Reference File Types
====================

The ``linearity`` step uses a LINEARITY reference file.

LINEARITY Reference File
-------------------------

:REFTYPE: LINEARITY
:Data model: `~jwst.datamodels.LinearityModel`

The LINEARITY reference file contains pixel-by-pixel polynomial correction
coefficients.

.. include:: linearity_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for LINEARITY
+++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in LINEARITY reference files,
because they are used as CRDS selectors
(see :ref:`linearity_selectors`):

=========  ==============================  ===========
Keyword    Data Model Name                 Instruments
=========  ==============================  ===========
DETECTOR   model.meta.instrument.detector  All
SUBARRAY   model.meta.subarray.name        All
FILTER     model.meta.instrument.filter    MIRI only
BAND       model.meta.instrument.band      MIRI only
=========  ==============================  ===========

Reference File Format
+++++++++++++++++++++
LINEARITY reference files are FITS format, with 2 IMAGE extensions
and 1 BINTABLE extension. The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  =======================  =========
EXTNAME  XTENSION  NAXIS  Dimensions               Data type
=======  ========  =====  =======================  =========
COEFFS   IMAGE       3    ncols x nrows x ncoeffs  float
DQ       IMAGE       2    ncols x nrows            integer
DQ_DEF   BINTABLE    2    TFIELDS = 4              N/A
=======  ========  =====  =======================  =========

Each plane of the COEFFS data cube contains the pixel-by-pixel coefficients for
the associated order of the polynomial. There can be any number of planes to
accommodate a polynomial of any order.

.. include:: ../includes/dq_def.rst
