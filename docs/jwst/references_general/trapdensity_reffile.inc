.. _trapdensity_reffile:

TRAPDENSITY Reference File
--------------------------

:REFTYPE: TRAPDENSITY
:Data model: `~jwst.datamodels.TrapDensityModel`

The TRAPDENSITY reference file contains a pixel-by-pixel map of the
trap density.

.. include:: ../references_general/trapdensity_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for TRAPDENSITY
++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in TRAPDENSITY reference files,
because they are used as CRDS selectors
(see :ref:`trapdensity_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
=========  ==============================

Reference File Format
+++++++++++++++++++++
TRAPDENSITY reference files are FITS format, with 2 IMAGE extensions
and 1 BINTABLE extension.
The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   int
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

.. include:: ../includes/dq_def.inc
