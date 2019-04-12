:orphan:

.. _ifufore_reffile:

IFUFORE Reference File (NIRSpec only)
-------------------------------------

:REFTYPE: IFUFORE
:Data model: `~jwst.datamodels.IFUFOREModel`

Reference Selection Keywords for IFUFORE
++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate IFUFORE references based on the following keywords.
IFUFORE is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ======================================
Instrument Keywords
========== ======================================
NIRSpec    INSTRUME, EXP_TYPE, DATE-OBS, TIME-OBS
========== ======================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
This file provides the parameters (Paraxial and distortions coefficients)
for the coordinate transforms from the MSA plane to the plane of the IFU slicer.

:model: Compound model, Polynomials

