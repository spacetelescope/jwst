.. _wavecorr_reffile:

WAVECORR Reference File
-----------------------

:REFTYPE: WAVECORR
:Data model: `~jwst.datamodels.WaveCorrModel`

The WAVECORR reference file contains pixel offset values as a function of
wavelength and source offset within a NIRSpec slit.
It is used when applying the NIRSpec wavelength zero-point correction to
fixed-slit (EXP_TYPE="NRS_FIXEDSLIT"), bright object TSO
(EXP_TYPE="NRS_BRIGHTOBJ"), and MSA/MOS spectra (EXP_TYPE="NRS_MSASPEC").
This is an optional correction that is turned on by default.
It can be turned off by specifying ``apply_wavecorr=False`` when running the step.

.. include:: ../references_general/wavecorr_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for WAVECORR
+++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in WAVECORR reference files,
because they are used as CRDS selectors
(see :ref:`wavecorr_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
+++++++++++++++++++++
WAVECORR reference files are in ASDF format, with the format and contents
specified by the `~jwst.datamodels.WaveCorrModel` data model schema.

