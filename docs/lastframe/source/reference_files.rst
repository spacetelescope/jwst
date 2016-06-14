Reference File
==============
The last frame correction step uses a LastFrame reference file.

CRDS Selection Criteria
-----------------------
Last frame reference files are selected on the basis of INSTRUME, DETECTOR, 
READPATT, and SUBARRAY values for the input science data set.

LastFrame Reference File Format
-------------------------------
The last frame correction reference files are FITS files with 3 
IMAGE extensions and 1 BINTABLE extension. The FITS primary data array 
is assumed to be empty. The  characteristics of the three image extension
are as follows:

=======  =====  =============  =========
EXTNAME  NAXIS  Dimensions     Data type
=======  =====  =============  =========
SCI      2      ncols x nrows  float
ERR      2      ncols x nrows  float
DQ       2      ncols x nrows  integer
=======  =====  =============  =========

The BINTABLE extension contains the bit assignments used in the DQ array.
It uses ``EXTNAME=DQ_DEF`` and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition

The SCI, ERR data and DQ  arrays are 2-D, with dimensions of ncols x nrows,
where ncols x nrows matches the dimensions of the raw detector
readout mode for which the last frame correction applies. 
