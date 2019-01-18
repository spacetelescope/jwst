:orphan:

.. _distortion_reffile:

DISTORTION Reference File
-------------------------

:REFTYPE: DISTORTION
:Data model: `~jwst.datamodels.DistortionModel`, `~jwst.datamodels.DistortionMRSModel`

Reference Selection Keywords for DISTORTION
+++++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate DISTORTION references based on the following keywords.
DISTORTION is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== =================================================================
Instrument Keywords
========== =================================================================
FGS        INSTRUME, DETECTOR, EXP_TYPE, DATE-OBS, TIME-OBS
MIRI       INSTRUME, DETECTOR, EXP_TYPE, CHANNEL, BAND, DATE-OBS, TIME-OBS
NIRCam     INSTRUME, DETECTOR, EXP_TYPE, CHANNEL, FILTER, DATE-OBS, TIME-OBS
NIRISS     INSTRUME, EXP_TYPE, PUPIL, DATE-OBS, TIME-OBS
========== =================================================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The distortion reference file contains a combination of astropy models,
representing the transform from detector to the telescope V2, V3 system.
The following convention was adopted:

- The output in the V2, V3 system is in units of arcsec.
- The input x and y are 0-based coordinates in the DMS system.
- The center of the first pixel is (0, 0), so the first pixel goes from -0.5 to 0.5.
- The origin of the transform is taken to be (0, 0).
  Note, that while a different origin can be used  for some transforms the relevant
  offset should first be prepended to the distortion transform to account for the change
  in origin of the coordinate frame.  For instance, MIRI takes input in (0, 0) - indexed
  detector pixel coordinates, but shifts these around prior to calling transforms that are
  defined with respect to science-frame pixels that omit reference pixels.


Internally the WCS pipeline works with 0-based coordinates.
When FITS header keywords are used, the 1 pixel offset in FITS coordinates is accounterd for
internally in the pipeline.

The model is a combination of polynomials.

:model: Transform from detector to an intermediate frame (instrument dependent).
