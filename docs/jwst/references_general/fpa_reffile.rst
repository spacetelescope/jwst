:orphan:

.. _fpa_reffile:

FPA Reference File (NIRSpec only)
---------------------------------

:REFTYPE: FPA
:Data model: `~jwst.datamodels.FPAModel`

Reference Selection Keywords for FPA
++++++++++++++++++++++++++++++++++++
CRDS selects appropriate FPA references based on the following keywords.
FPA is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ======================================
Instrument Keywords
========== ======================================
NIRSpec    INSTRUME, EXP_TYPE, DATE-OBS, TIME-OBS
========== ======================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The FPA reference file stores information on the metrology of the Focal Plane Assembly (FPA),
which consists of two Sensor Chip Arrays (SCA), named NRS1 and NRS2.

The reference file contains two fields : “nrs1_model” and “nrs2_model”.
Each of them stores the transform (shift and rotation) to transform positions
from the FPA to the respective SCA. The output units are in pixels.

:nrs1_model: Transform for the NRS1 detector.
:nrs2_model: Transform for the NRS2 detector.

