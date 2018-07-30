Reference File Types
====================
The Data Quality Initialization step uses a MASK reference file.

CRDS Selection Criteria
-----------------------
MASK reference files are currently selected based only on the value of
DETECTOR in the input science data set. There is one MASK reference file for
each JWST instrument detector.

MASK Reference File Format
--------------------------
The MASK reference file is a FITS file with a primary HDU, 1 IMAGE extension
HDU and 1 BINTABLE extension. The primary data array is assumed to be empty.
The MASK data are stored in the first IMAGE extension, which shall have
EXTNAME='DQ'. The data array in this extension has integer data type and is
2-D, with dimensions equal to the number of columns and rows in a full frame
raw readout for the given detector, including reference pixels. Note that
this does not include the reference output for MIRI detectors.

The BINTABLE extension contains the bit assignments used in the DQ array.
It uses ``EXTNAME=DQ_DEF`` and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition


