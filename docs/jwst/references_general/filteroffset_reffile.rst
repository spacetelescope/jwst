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

.. include:: ../includes/standard_keywords.rst

Reference File Format
+++++++++++++++++++++
The filter offset reference file is an ASDF file that contains a dictionary of row and column offsets for the MIRI imaging dataset. The filter offset reference file contains a dictionary in the tree that is indexed by the instrument filter. Each filter points to two fields - row_offset and column_offset. The format is

:miri_filter_name:
    :column_offset: Offset in x (in arcmin)
    :row_offset: Offset in y (in arcmin)

