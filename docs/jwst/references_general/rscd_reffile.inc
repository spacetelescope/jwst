.. _rscd_reffile:

RSCD Reference File
-------------------

:REFTYPE: RSCD
:Data model: `~jwst.datamodels.RSCDModel`

The RSCD reference file contains the number of groups to flag as 'DO_NOT_USE' based on
readout mode and subarray size. Integrations two and higher are flagged; the first integration
groups values are not flagged. 

.. include:: ../references_general/rscd_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for RSCD
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in RSCD reference files,
because they are used as CRDS selectors
(see :ref:`rscd_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
=========  ==============================

Reference File Format
+++++++++++++++++++++
RSCD reference files are in FITS format. The number of groups to set to 'DO_NOT_USE'
is contained in the first BINTABLE extension.
The FITS primary HDU does not contain a data array.
The BINTABLE extension uses the identifier EXTNAME = "RSCD_GROUP_SKIP" and
the characteristics of the table columns are as follows:

===========  ==========  ==================================================
Column name  Data type   Notes
===========  ==========  ==================================================
subarray     char\*13    FULL or subarray name
readpatt     char\*8     SLOWR1, FASTR1, SLOW, or FAST
group_skip   int         number of initial groups in an integration to flag
===========  ==========  ==================================================

The entries in the first two columns of the table are used as selection criteria, matching
the exposure properties of the data. The last column contains the
number of initial groups in an integration to flag as 'DO_NOT_USE' for integrations 2 and higher.
These initial groups are affected by the RSCD effect.  The groups in the first integration are not flagged. 

