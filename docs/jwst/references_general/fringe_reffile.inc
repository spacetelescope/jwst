.. _fringe_reffile:

FRINGE Reference File
---------------------

:REFTYPE: FRINGE
:Data model: `~jwst.datamodels.FringeModel`

The FRINGE reference file contains pixel-by-pixel fringing correction
values.

.. include:: ../references_general/fringe_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for FRINGE
+++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in FRINGE reference files,
because they are used as CRDS selectors
(see :ref:`fringe_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
BAND       model.meta.instrument.band
=========  ==============================

Reference File Format
+++++++++++++++++++++
FRINGE reference files are FITS format, with 3 IMAGE extensions
and 1 BINTABLE extension. The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
ERR      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   integer
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

The values in the ``SCI`` array give the correction values to be applied to
the science data.
Because MIRI MRS exposures are always full-frame, the image dimensions
should be 1032 x 1024.

.. include:: ../includes/dq_def.inc

