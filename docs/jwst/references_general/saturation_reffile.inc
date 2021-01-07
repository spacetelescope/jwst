.. _saturation_reffile:

SATURATION Reference File
-------------------------

:REFTYPE: SATURATION
:Data model: `jwst.datamodels.SaturationModel`

The SATURATION reference file contains pixel-by-pixel saturation threshold
values.

.. include:: ../references_general/saturation_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for SATURATION
+++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in SATURATION reference files,
because they are used as CRDS selectors
(see :ref:`saturation_selectors`):

=========  ==============================  ==========
Keyword    Data Model Name                 Instrument
=========  ==============================  ==========
DETECTOR   model.meta.instrument.detector  All
SUBARRAY   model.meta.subarray.name        NIRSpec
=========  ==============================  ==========

Reference File Format
+++++++++++++++++++++
SATURATION reference files are FITS format, with 2 IMAGE extensions
and 1 BINTABLE extension. The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   integer
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

The values in the ``SCI`` array give the saturation threshold in units of DN
for each pixel.

.. include:: ../includes/dq_def.inc
