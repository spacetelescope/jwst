:orphan:

.. _filteroffset_reffile:


FILTEROFFSET Reference File (MIRI only)
---------------------------------------

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
========== ================================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The filteroffset reference file is an ASDF file that contains a list
called ``filters``. Every item in the list contains three fields -
``row_offset``, ``column_offset`` and ``name`` - the name of the filter they are valid for. The offsets, in pixels, are applied in the image science frame.

:filters:
    :column_offset: Offset in x (in pix)
    :row_offset: Offset in y (in pix)
    :name: Filter name
