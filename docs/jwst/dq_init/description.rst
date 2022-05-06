Description
============

:Class: `jwst.dq_init.DQInitStep`
:Alias: dq_init

The Data Quality (DQ) initialization step in the calibration pipeline
populates the DQ mask for the input dataset. Flag values from the
appropriate static mask ("MASK") reference file in CRDS are copied into the
"PIXELDQ" array of the input dataset, because it is assumed that flags in the
mask reference file pertain to problem conditions that affect all groups and
integrations for a given pixel.

The actual process consists of the following steps:

 - Determine what MASK reference file to use via the interface to the bestref
   utility in CRDS.

 - If the "PIXELDQ" or "GROUPDQ" arrays of the input dataset do not already exist,
   which is sometimes the case for raw input products, create these arrays in
   the input data model and initialize them to zero. The "PIXELDQ" array will be
   2D, with the same number of rows and columns as the input science data.
   The "GROUPDQ" array will be 4D with the same dimensions (nints, ngroups,
   nrows, ncols) as the input science data array.

 - Check to see if the input science data is in subarray mode. If so, extract a
   matching subarray from the full-frame MASK reference file.

 - Propagate the DQ flags from the reference file DQ array to the science data "PIXELDQ"
   array using numpy's ``bitwise_or`` function.

Note that when applying the ``dq_init`` step to FGS guide star data, as is done in
the :ref:`calwebb_guider <calwebb_guider>` pipeline, the flags from the MASK reference
file are propagated into the guide star dataset "DQ" array, instead of the "PIXELDQ" array.
The step identifies guide star data based on the following exposure type (EXP_TYPE keyword) values:
FGS_ID-IMAGE, FGS_ID-STACK, FGS_ACQ1, FGS_ACQ2, FGS_TRACK, and FGS_FINEGUIDE.

NIRSpec IRS2
------------
No special handling is required for NIRSpec exposures taken using the IRS2
readout pattern, because matching IRS2 MASK reference files are supplied in CRDS.
