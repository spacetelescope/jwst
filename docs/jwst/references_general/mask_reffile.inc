.. _mask_reffile:

MASK Reference File
-------------------

:REFTYPE: MASK
:Data model: `~jwst.datamodels.MaskModel`

The MASK reference file contains pixel-by-pixel DQ flag values that indicate
problem conditions.

.. include:: ../references_general/mask_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for MASK
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in MASK reference files,
because they are used as CRDS selectors
(see :ref:`mask_selectors`):

=========  ==============================  ============
Keyword    Data Model Name                 Instruments
=========  ==============================  ============
DETECTOR   model.meta.instrument.detector  All
SUBARRAY   model.meta.subarray.name        All
EXP_TYPE   model.meta.exposure.type        FGS only
READPATT   model.meta.exposure.readpatt    NIRSpec only
=========  ==============================  ============

Reference File Format
+++++++++++++++++++++
MASK reference files are FITS format, with one IMAGE extension
and 1 BINTABLE extension. The FITS primary HDU does not contain a
data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
DQ       IMAGE       2    ncols x nrows   integer
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

The values in the ``DQ`` array give the per-pixel flag conditions that are
to be propagated into the science exposure's ``PIXELDQ`` array.
The dimensions of the ``DQ`` array should be equal to the number of columns
and rows in a full-frame readout of a given detector, including reference
pixels.
Note that this does not include the reference *output* for MIRI detectors.

.. include:: ../includes/dq_def.inc
