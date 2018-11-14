Description
============
The Data Quality (DQ) initialization step in the calibration pipeline
populates the DQ mask for the input dataset. Flags from the
appropriate static mask reference file in CRDS are copied into the
``PIXELDQ`` array of the input dataset, because it is assumed that flags in the
mask reference file pertain to problem conditions that are group- and
integration-independent.

The actual process consists of the following steps:

 - Determine what MASK reference file to use via the interface to the bestref
   utility in CRDS.

 - If the ``PIXELDQ`` or ``GROUPDQ`` arrays of the input dataset do not already exist,
   which is usually the case for raw input products, create these arrays in
   the input data model and initialize them to zero. The ``PIXELDQ`` array will be
   2-D, with the same number of rows and columns as the input science data.
   The ``GROUPDQ`` array will be 4-D with the same dimensions (nints, ngroups,
   nrows, ncols) as the input science data array.

 - Check to see if the input science data is in subarray mode. If so, extract a
   matching subarray from the full-frame MASK reference file.

 - Copy the DQ flags from the reference file mask to the science data ``PIXELDQ``
   array using numpy's ``bitwise_or`` function.
