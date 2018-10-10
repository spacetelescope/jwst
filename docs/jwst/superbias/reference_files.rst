Reference File Types
====================
The superbias subtraction step uses a SUPERBIAS reference file.

.. include:: ../includes/standard_keywords.rst

.. include:: superbias_selection.rst

SUPERBIAS Reference File Format
-------------------------------
Superbias reference files are FITS files with 3 IMAGE extensions and 1 BINTABLE
extension. The FITS primary data array is assumed to be empty. The 
characteristics of the three image extension are as follows:

=======  =====  =============  =========
EXTNAME  NAXIS  Dimensions     Data type
=======  =====  =============  =========
SCI      2      ncols x nrows  float
ERR      2      ncols x nrows  float
DQ       2      ncols x nrows  integer
=======  =====  =============  =========

.. include:: ../includes/dq_def.rst
