.. _bkg_reffile:

Background reference file
-------------------------

:REFTYPE: BKG

The WFSS background reference file contains a "master" image of the
dispersed background produced by a particular filter+grism combination.

.. include:: ../references_general/bkg_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for BKG
++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in BKG reference files,
because they are used as CRDS selectors.
(see :ref:`bkg_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
EXP_TYPE   model.meta.exposure.type
FILTER     model.meta.instrument.filter
PUPIL      model.meta.instrument.pupil
=========  ==============================

Reference File Format
+++++++++++++++++++++
BKG reference files are in FITS format. The content of the background references differ
based on exposure type:

:Data model: `~jwst.datamodels.WfssBkgModel`

Wide-field slitless spectroscopy (WFSS) background files contain 3 IMAGE extensions and
1 BINTABLE extension. The FITS primary HDU does not contain a data array.
The characteristics of the FITS extensions are as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
ERR      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   integer
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

.. include:: ../includes/dq_def.inc

:Data model: `~jwst.datamodels.SossBkgModel`

Background reference files supporting NIRISS SOSS exposures contain three
image extensions, each with three dimensions. The extensions contain templates
of reference backgrounds. The FITS primary HDU does not contain a data array.
The characteristics of the FITS extensions are as follows:

=======  ========  =====  ==========================   =========
EXTNAME  XTENSION  NAXIS  Dimensions                   Data type
=======  ========  =====  ==========================   =========
SCI      IMAGE       3    ncols x nrows x ntemplates   float
ERR      IMAGE       3    ncols x nrows x ntemplates   float
DQ       IMAGE       3    ncols x nrows x ntemplates   integer
=======  ========  =====  ==========================   =========
