Reference File Types
====================
The saturation step uses a SATURATION reference file.

.. include:: ../includes/standard_keywords.rst

.. include:: saturation_selection.rst
             

SATURATION Reference File Format
--------------------------------
Saturation reference files are FITS format with
with 2 IMAGE extensions: ``SCI`` and ``DQ``, which are both 2-D integer arrays,
and 1 BINTABLE extension.

The values in the SCI array give the saturation threshold in units of DN for
each pixel. The saturation reference file also contains a ``DQ_DEF`` table
extension, which lists the bit assignments for the flag conditions used in
the DQ array.

.. include:: ../includes/dq_def.rst

