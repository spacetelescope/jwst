:orphan:

.. _ifupost_reffile:
  
IFUPOST Reference File (NIRSpec only)
-------------------------------------

:REFTYPE: IFUPOST
:Data model: `~jwst.datamodels.IFUPostModel`

Reference Selection Keywords for IFUPOST
++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate IFUPOST references based on the following keywords.
IFUPOST is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ==============================================
Instrument Keywords
========== ==============================================
NIRSpec    INSTRUME, EXP_TYPE, OPMODE, DATE-OBS, TIME-OBS
========== ==============================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The IFUPOST reference file provides the parameters (Paraxial and distortions coefficients)
for the coordinate transforms from the slicer plane to the MSA plane (out),
that is the plane of the IFU virtual slits.

The reference file contains models made up based on an offset and a polynomial.
There is a model for each of the slits and is indexed by the slit number.
The models is used as part of the conversion from the GWA to slit.

:slice_<slice_number>:
    :model: Polynomial and rotation models.
