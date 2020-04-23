.. _abvega_offset_reffile:

ABVEGA_OFFSET Reference File
----------------------------

:REFTYPE: ABVEGA_OFFSET
:Data model: `~jwst.datamodels.ABVegaOffsetModel`

The ABVEGA_OFFSET reference file contains data necessary for converting
from AB to Vega magnitudes.

.. include:: ../references_general/abvega_offset_selection.inc

.. include:: ../includes/standard_keywords.inc

ABVEGA_OFFSET Reference File Format
+++++++++++++++++++++++++++++++++++

ABVEGA_OFFSET reference files are in ASDF format.  The ABVEGA_OFFSET
reference file contains tabular data in a key called
``abvega_offset``.  The content of the table varies for different
instrument modes, as shown in the tables below.

+------------+---------------+-----------+------------+----------+
| Instrument | Column name   | Data type | Dimensions | Units    |
+============+===============+===========+============+==========+
| FGS        | detector      | string    | 7          | N/A      |
+            +---------------+-----------+------------+----------+
|            | abvega_offset | float     | scalar     | unitless |
+------------+---------------+-----------+------------+----------+

+------------+---------------+-----------+------------+----------+
| Instrument | Column name   | Data type | Dimensions | Units    |
+============+===============+===========+============+==========+
| MIRI       | filter        | string    | 12         | N/A      |
+            +---------------+-----------+------------+----------+
|            | abvega_offset | float     | scalar     | unitless |
+------------+---------------+-----------+------------+----------+

+------------+---------------+-----------+------------+----------+
| Instrument | Column name   | Data type | Dimensions | Units    |
+============+===============+===========+============+==========+
| NIRCam     | filter        | string    | 12         | N/A      |
+ or         +---------------+-----------+------------+----------+
| NIRISS     | pupil         | string    | 15         | N/A      |
+            +---------------+-----------+------------+----------+
|            | abvega_offset | float     | scalar     | unitless |
+------------+---------------+-----------+------------+----------+


Row Selection
^^^^^^^^^^^^^

A row of data within the reference table is selected by the pipeline
step based on the optical elements in use for the exposure. The
selection attributes are always contained in the first few columns of
the table.  The last column contains the data needed to convert from
AB to Vega magnitudes.  The row selection criteria for each
instrument/mode are:

* FGS:
   - Detector
* MIRI:
   - Filter
* NIRCam:
   - Filter and Pupil
* NIRISS:
   - Filter and Pupil
