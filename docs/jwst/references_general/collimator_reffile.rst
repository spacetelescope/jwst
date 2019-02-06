:orphan:

.. _collimator_reffile:

COLLIMATOR Reference File (NIRSpec only)
----------------------------------------

:REFTYPE: COLLIMATOR
:Data model: `~jwst.datamodels.CollimatorModel`

Reference Selection Keywords for COLLIMATOR
+++++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate COLLIMATOR references based on the following keywords.
COLLIMATOR is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ======================================
Instrument Keywords
========== ======================================
NIRSpec    INSTRUME, EXP_TYPE, DATE-OBS, TIME-OBS
========== ======================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The collimator reference file contains an astropy compound model made up of polynomial models,
rotations, and translations. The forward direction is from the GWA to the MSA.

:model: Transform through the COLLIMATOR.
