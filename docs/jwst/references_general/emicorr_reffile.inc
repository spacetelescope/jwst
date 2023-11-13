.. _emicorr_reffile:

EMICORR Reference File
----------------------

:REFTYPE: EMI
:Data model: `~jwst.datamodels.EmiModel`

The EMI reference file contains data necessary for removing
contaminating MIRI EMI frequencies.

.. include:: ../references_general/emicorr_selection.inc

.. include:: ../includes/standard_keywords.inc

EMICORR Reference File Format
+++++++++++++++++++++++++++++

MIRI EMIVORR reference files are in ASDF format.  The EMICORR
reference file contains data in the following keyes:
``Hz390``, ``Hz218a``, ``Hz218b``, ``Hz218c``, ``Hz164`` and
``Hz10``.  The content of each entry is a 1-D array with the
phase amplitude values corresponding to each frequency.

+------------+---------------+-----------+------------+----------+
| Instrument | Key name      | Data type | Dimensions | Units    |
+============+===============+===========+============+==========+
| MIRI       | Hz390         | float     | scalar     | Hz       |
+            +---------------+-----------+------------+----------+
|            | Hz218a        | float     | scalar     | Hz       |
+            +---------------+-----------+------------+----------+
|            | Hz218b        | float     | scalar     | Hz       |
+            +---------------+-----------+------------+----------+
|            | Hz218c        | float     | scalar     | Hz       |
+            +---------------+-----------+------------+----------+
|            | Hz164         | float     | scalar     | Hz       |
+            +---------------+-----------+------------+----------+
|            | Hz10          | float     | scalar     | Hz       |
+------------+---------------+-----------+------------+----------+


Frequency Selection
^^^^^^^^^^^^^^^^^^^

The frequency to be corrected for will be selected according to the
dictionary ``subarray_cases``, key ``freqs``, located in the
``emicorr.py`` script within the ``emicorr` directory under ``jwst``.
