:orphan:

.. _fore_reffile:

FORE Reference File (NIRSpec only)
----------------------------------

:REFTYPE: FORE
:Data model: `~jwst.datamodels.FOREModel`

Reference Selection Keywords for FORE
+++++++++++++++++++++++++++++++++++++
CRDS selects appropriate FORE references based on the following keywords.
FORE is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ==============================================
Instrument Keywords
========== ==============================================
NIRSpec    INSTRUME, EXP_TYPE, FILTER, DATE-OBS, TIME-OBS
========== ==============================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The FORE reference file stores the transform through the Filter Wheel Assembly (FWA).
It has two fields - “filter” and “model”. The transform through the FWA is chromatic.
It is represented as a Polynomial of two variables whose coefficients are wavelength dependent.
The compound model takes three inputs - x, y positions and wavelength.

:filter: Filter name.
:model: Transform through the Filter Wheel Assembly (FWA).

