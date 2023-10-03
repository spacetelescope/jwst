:orphan:

.. _filteroffset_reffile:


FILTEROFFSET Reference File
---------------------------

:REFTYPE: FILTEROFFSET
:Data model: `~jwst.datamodels.FilteroffsetModel`

Reference Selection Keywords for FILTEROFFSET
+++++++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate FILTEROFFSET references based on the following keywords.
FILTEROFFSET is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ================================================
Instrument Keywords
========== ================================================
MIRI       INSTRUME, DETECTOR, EXP_TYPE, DATE-OBS, TIME-OBS
NIRCam     INSTRUME, CHANNEL, MODULE, DATE-OBS, TIME-OBS
NIRISS     INSTRUME, EXP_TYPE, DATE-OBS, TIME-OBS
========== ================================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The filteroffset reference file is an ASDF file that contains a list
called ``filters``. Every item in the list contains one or more entries that
are used as selectors, as well as the column and row offset values to be applied.
For the MIRI instrument, there is one selector "filter", which is the name of
the filter to which the offsets apply. For NIRCam and NIRISS, the selectors are "filter" and
"pupil".

The offsets, in units of pixels, are *added* to positions in images that use
the reference filter/pupil, in order to align images to the reference filter/pupil
frame.

:filters:
    :filter: Filter name
    :pupil: Pupil name (NIRCam only)
    :column_offset: Offset in x (in pixels)
    :row_offset: Offset in y (in pixels)
