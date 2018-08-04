Reference File Types
====================
The saturation step uses a SATURATION reference file.

CRDS Selection Criteria
-----------------------
Saturation reference files are selected on the basis of INSTRUME, DETECTOR, and 
SUBARRAY values from the input science data set.

SATURATION Reference File Format
--------------------------------
Saturation reference files are FITS format with
with 2 IMAGE extensions: ``SCI`` and ``DQ``, which are both 2-D integer arrays,
and 1 BINTABLE extension.

The values in the SCI array give the saturation threshold in units of DN for
each pixel. The saturation reference file also contains a ``DQ_DEF`` table
extension, which lists the bit assignments for the flag conditions used in
the DQ array.

The BINTABLE extension uses ``EXTNAME=DQ_DEF`` and contains the bit assignments
of the conditions flagged in the DQ array, and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition

