Reference File
==============
The dark current step uses a DARK reference file.

CRDS Selection Criteria
-----------------------
Dark reference files are selected on the basis of INSTRUME, DETECTOR, and 
SUBARRAY values for the input science data set. For MIRI exposures, the
value of READPATT is used as an additional selection criterion.

DARK Reference File Format
--------------------------
Dark reference files are FITS files with 3 IMAGE extensions and 1 BINTABLE
extension. The FITS primary data array is assumed to be empty. The 
characteristics of the three image extension are as follows:

=======  =====  =======================  =========
EXTNAME  NAXIS  Dimensions               Data type
=======  =====  =======================  =========
SCI      3      ncols x nrows x ngroups  float
ERR      3      ncols x nrows x ngroups  float
DQ       2      ncols x nrows            integer
=======  =====  =======================  =========

The BINTABLE extension contains the bit assignments used in the DQ array.
It uses ``EXTNAME=DQ_DEF`` and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition
