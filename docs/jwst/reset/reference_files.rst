Reference File Types
====================
The reset correction step uses a RESET reference file.

CRDS Selection Criteria
-----------------------
Reset reference files are selected on the basis of INSTRUME, DETECTOR,
READPATT and SUBARRAY values for the input science data set.

RESET Reference File Format
---------------------------
The reset reference files are FITS files with 3 IMAGE extensions and 1 BINTABLE
extension. The FITS primary data array is assumed to be empty. The 
characteristics of the three image extension are as follows:

=======  =====  ============================== =========
EXTNAME  NAXIS  Dimensions                     Data type
=======  =====  ============================== =========
SCI      4      ncols x nrows x ngroups x nint float
ERR      4      ncols x nrows x ngroups x nint float
DQ       2      ncols x nrows                  integer
=======  =====  ============================== =========

The BINTABLE extension contains the bit assignments used in the DQ array.
It uses ``EXTNAME=DQ_DEF`` and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition


The SCI and ERR data arrays are 4-D, with dimensions of ncols x nrows x 
ngroups X nints, where ncols x nrows matches the dimensions of the raw detector
readout mode for which the reset applies. The reference file contains the 
number of NGroups planes required for the correction to be zero on
the last plane Ngroups plane.  The correction for the first few
integrations varies and eventually settles down to a constant correction
independent of integration number.  
