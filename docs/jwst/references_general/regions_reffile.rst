:orphan:

.. _regions_reffile:
  
REGIONS Reference File (MIRI only)
----------------------------------

:REFTYPE: REGIONS
:Data model: `~jwst.datamodels.RegionsModel`

Reference Selection Keywords for REGIONS
++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate REGIONS references based on the following keywords.
REGIONS is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ===============================================================
Instrument Keywords
========== ===============================================================
MIRI       INSTRUME, DETECTOR, EXP_TYPE, CHANNEL, BAND, DATE-OBS, TIME-OBS
========== ===============================================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The IFU takes a region reference file that defines the pixels on the detector that map to the sky and which
IFU slice they correspond to.  The reference file is a three-dimensional array
measuring 1032 x 1024 x 9.  Each of the nine planes within the array represents
a two-dimensional detector image for which each pixel has a value indicating the slice to which it belongs.
The nine planes each correspond to different throughput thresholds ranging from 10% - 90%; these differ slightly
in the effective along-slice size.

:channel: The MIRI channels in the observation, e.g. "12".
:band: The band for the observation (one of "LONG", "MEDIUM", "SHORT").
:regions: An array with the size of the MIRI MRS image where pixel values map to the MRS slice number. 0 indicates a pixel is not within any slice.
