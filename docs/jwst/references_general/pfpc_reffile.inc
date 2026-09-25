.. _pfpc_reffile:

PFPC Reference File
---------------------------

:REFTYPE: PFPC

The PFPC reference file contains data necessary for applying a point fixed
pattern correction to dithered spectral data.

.. include:: ../includes/standard_keywords.inc


PFPC Reference File Format
++++++++++++++++++++++++++++++++++

PFPC reference files are in FITS format with one table extension with ``EXTNAME = 'PFPC'``.
The primary HDU does not contain a data array.  The contents of the table vary for
different instrument modes, as shown in the tables below.  All tables must include
``wavelength`` and ``correction`` columns.  All other columns should be FITS keywords,
used for matching the correct row to the input data.

.. todo: add reference to model when it is merged

:Data model: stdatamodels.jwst.datamodels.MirMrsPFPCModel

+------------+-------+----------------+-----------+------------+-------+
| Instrument | Mode  | Column name    | Data type | Dimensions | Units |
+============+=======+================+===========+============+=======+
| MIRI       | MRS   | channel        | string    | scalar     | N/A   |
+            +       +----------------+-----------+------------+-------+
|            |       | band           | string    | scalar     | N/A   |
+            +       +----------------+-----------+------------+-------+
|            |       | patttype       | string    | scalar     | N/A   |
+            +       +----------------+-----------+------------+-------+
|            |       | dithdirc       | string    | scalar     | N/A   |
+            +       +----------------+-----------+------------+-------+
|            |       | dithopfr       | string    | scalar     | N/A   |
+            +       +----------------+-----------+------------+-------+
|            |       | mrsprchn       | string    | scalar     | N/A   |
+            +       +----------------+-----------+------------+-------+
|            |       | pattnum        | integer   | scalar     | N/A   |
+            +       +----------------+-----------+------------+-------+
|            |       | r_flat         | string    | scalar     | N/A   |
+            +       +----------------+-----------+------------+-------+
|            |       | r_photom       | string    | scalar     | N/A   |
+            +       +----------------+-----------+------------+-------+
|            |       | wavelength     | float     | 1D vector  | um    |
+            +       +----------------+-----------+------------+-------+
|            |       | correction     | float     | 1D vector  | N/A   |
+------------+-------+----------------+-----------+------------+-------+
