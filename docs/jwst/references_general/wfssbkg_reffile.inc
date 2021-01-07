.. _wfssbkg_reffile:

WFSS Background reference file
------------------------------

:REFTYPE: WFSSBKG
:Data model: `~jwst.datamodels.WfssBkgModel`

The WFSS background reference file contains a "master" image of the
dispersed background produced by a particular filter+grism combination.

.. include:: ../references_general/wfssbkg_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for WFSSBKG
++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in WFSSBKG reference files,
because they are used as CRDS selectors
(see :ref:`wfssbkg_selectors`):

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
WFSSBKG reference files are FITS files with 3 IMAGE extensions and
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
