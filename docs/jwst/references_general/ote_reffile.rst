:orphan:

.. _ote_reffile:
  
OTE Reference File (NIRSpec only)
---------------------------------

:REFTYPE: OTE
:Data model: `~jwst.datamodels.OTEModel`

Reference Selection Keywords for OTE
++++++++++++++++++++++++++++++++++++
CRDS selects appropriate OTE references based on the following keywords.
OTE is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ======================================
Instrument Keywords
========== ======================================
NIRSpec    INSTRUME, EXP_TYPE, DATE-OBS, TIME-OBS
========== ======================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The OTE reference file contains a combination of astropy models - polynomial, shift, rotation and scaling.

:model: Transform through the Optical Telescope Element (OTE), from the FWA to XAN, YAN telescope frame. The
        output units are in arcsec.
