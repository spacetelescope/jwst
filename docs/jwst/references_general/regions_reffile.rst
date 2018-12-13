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

.. include:: ../includes/standard_keywords.rst

Reference File Format
+++++++++++++++++++++
The IFU takes a region reference file that defines the region over which the WCS is valid.
The reference file should define a polygon and may consist of a set of X,Y coordinates that define the polygon.

:channel: The MIRI channels in the observation, e.g. "12".
:band: The band for the observation (one of "LONG", "MEDIUM", "SHORT").
:regions: An array with the size of the MIRI MRS image where pixel values map to the MRS slice number. 0 indicates a pixel is not within any slice.
