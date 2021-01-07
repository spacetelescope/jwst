.. _flat_reffile:

FLAT Reference File
-------------------

:REFTYPE: FLAT
:Data model: `~jwst.datamodels.FlatModel`

The FLAT reference file contains pixel-by-pixel detector response values.
It is used for all instrument modes except the NIRSpec spectroscopic modes.

.. include:: ../references_general/flat_selection.inc

.. include:: ../includes/standard_keywords.inc

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

The ``DQ_DEF`` table extension lists the bit assignments for the flag conditions
used in the DQ array.

.. include:: ../includes/dq_def.inc

For application to imaging data, the FITS file contains a single set of SCI,
ERR, DQ, and DQ_DEF extensions.  Image dimensions should be 2048x2048 for the
NIR detectors and 1032x1024 for MIRI (i.e. they include reference pixels),
unless data were taken in subarray mode.

