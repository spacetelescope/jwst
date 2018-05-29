Reference File Types
====================

The linearity correction step uses a LINEARITY reference file.

CRDS Selection Criteria
-----------------------
Linearity reference files are selected by INSTRUME and DETECTOR.

Reference File Format
---------------------
Linearity reference files are FITS format with 2 IMAGE extensions and 1
BINTABLE extension. The primary data array is assumed to be empty. The 2
IMAGE extensions have the following characteristics:

=======  =====  =======================  =========
EXTNAME  NAXIS  Dimensions               Data type
=======  =====  =======================  =========
COEFFS   3      ncols x nrows x ncoeffs  float
DQ       2      ncols x nrows            integer
=======  =====  =======================  =========

Each plane of the COEFFS data cube contains the pixel-by-pixel coefficients for
the associated order of the polynomial. There can be any number of planes to
accommodate a polynomial of any order.

The BINTABLE extension uses ``EXTNAME=DQ_DEF`` and contains the bit assignments
of the conditions flagged in the DQ array.

