:orphan:

.. _camera_reffile:

CAMERA Reference File (NIRSpec only)
------------------------------------

:REFTYPE: CAMERA
:Data model: `~jwst.datamodels.CameraModel`

Reference Selection Keywords for CAMERA
+++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate CAMERA references based on the following keywords.
CAMERA is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ======================================
Instrument Keywords
========== ======================================
NIRSpec    INSTRUME, EXP_TYPE, DATE-OBS, TIME-OBS
========== ======================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The camera reference file contains an astropy compound model made up of polynomial models, rotations, and translations. The forward direction is from the FPA to the GWA.

:model: Transform through the CAMERA.

