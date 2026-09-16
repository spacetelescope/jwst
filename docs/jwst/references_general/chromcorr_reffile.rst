:orphan:

.. _chromcorr_reffile:

CHROMCORR Reference File (NIRSpec only)
----------------------------------------

:REFTYPE: CHROMCORR
:Data model: `~stdatamodels.jwst.datamodels.ChromCorrModel`

Reference Selection Keywords for CHROMCORR
++++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate CHROMCORR references based on the following keywords.
CHROMCORR is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ======================================
Instrument Keywords
========== ======================================
NIRSpec    INSTRUME, EXP_TYPE, FILTER
========== ======================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The chromatic correction reference file contains an astropy compound model made up of
polynomial models, translations, and other transformations.

:model: Chromaticity correction transform.
