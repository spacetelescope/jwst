Description
============
The Data Quality Initialization step in the calibration pipeline
populates the Data Quality mask for the input dataset. DQ flags from the
appropriate static mask reference file in CRDS are copied into the 
PIXELDQ object of the input dataset, because it is assumed that flags in the
mask reference file pertain to problem conditions that are group- and 
integration-independent.

The actual process consists of the following steps:

 - Determine what mask reference file to use via the interface to the bestref
   utility in CRDS.

 - If the PIXELDQ and GROUPDQ objects of the input dataset do not already exist,
   which is the case for raw Level-1b input products, create these objects in
   the input data model and initialize them to zero. The PIXELDQ array will be 
   2-D, with the same number of rows and columns as the input science data.
   The GROUPDQ array will be 4-D with the same dimensions (nints, ngroups,
   nrows, ncols) as the input science data array.

 - Check to see if the input science data is in subarray mode. If so, extract a
   matching subarray from the full frame mask reference file.

 - Copy the DQ flags from the reference file mask to the science data PIXELDQ
   array using numpy's bitwise_or function.

