.. _superbias_reffile:

SUPERBIAS Reference File
------------------------

:REFTYPE: SUPERBIAS
:Data model: `~jwst.datamodels.SuperBiasModel`

The SUPERBIAS reference file contains a 2-D image of the detector bias
("zeroth" read) structure.

.. include:: ../references_general/superbias_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for SUPERBIAS
++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in SUPERBIAS reference files,
because they are used as CRDS selectors
(see :ref:`superbias_selectors`):

=========  ==============================  ============================
Keyword    Data Model Name                 Instruments
=========  ==============================  ============================
DETECTOR   model.meta.instrument.detector  FGS, NIRCam, NIRISS, NIRSpec
READPATT   model.meta.exposure.readpatt    FGS, NIRCam, NIRISS, NIRSpec
SUBARRAY   model.meta.subarray.name        FGS, NIRCam, NIRISS, NIRSpec
SUBSTRT1   model.meta.subarray.xstart      NIRSpec only
SUBSTRT2   model.meta.subarray.ystart      NIRSpec only
SUBSIZE1   model.meta.subarray.xsize       NIRSpec only
SUBSIZE2   model.meta.subarray.ysize       NIRSpec only
=========  ==============================  ============================

Reference File Format
+++++++++++++++++++++
SUPERBIAS reference files are FITS format, with 3 IMAGE extensions
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

The ``SCI`` array contains the super-bias image of the detector. The
``ERR`` array contains uncertainties in the super-bias values and the
``DQ`` array contains data quality flags associated with the super-bias
image.

.. include:: ../includes/dq_def.inc

