Reference File Types
====================

The fringe correction step uses a FRINGE reference file, which has the same
format as the FLAT reference file.  This correction is applied only to MIRI 
MRS (IFU) mode exposures, which are always single full-frame 2-D images.

CRDS Selection Criteria
-----------------------
Fringe reference files are selected by DETECTOR and GRATNG14.

Reference File Format
---------------------
Fringe reference files are FITS format with 3 IMAGE extensions and 1
BINTABLE extension. The primary data array is assumed to be empty. The 3
IMAGE extensions have the following characteristics:

=======  =====  =============  =========
EXTNAME  NAXIS  Dimensions     Data type
=======  =====  =============  =========
SCI      2      ncols x nrows  float
ERR      2      ncols x nrows  float
DQ       2      ncols x nrows  integer
=======  =====  =============  =========

Image dimensions should be 1032 x 1024.

The BINTABLE extension uses ``EXTNAME=DQ_DEF`` and contains the bit assignments
of the conditions flagged in the DQ array.

