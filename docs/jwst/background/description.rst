
Description
===========
The background subtraction step in the calibration pipeline performs
image-from-image subtraction in order to accomplish subtraction of background
signal. The step takes as input one target exposure, to which the
subtraction will be applied, and a list of one or more background exposures.

If more than one background exposure is provided, they will be averaged
together before being applied to the target exposure.
The average background image is produced as follows:

 - The SCI arrays of all background exposures are averaged
 - The ERR arrays of all background exposures are summed in quadrature and
   then converted to an uncertainty in the mean
 - The DQ arrays of all background exposures are combined using a bitwise-OR
   operation

The average background exposure is then subtracted from the target exposure.
The subtraction consists of the following operations:

 - The SCI array of the average background is subtracted from the SCI
   array of the target

 - The ERR array of the target is not operated upon until full error
   propagation is implemented in the entire pipeline

 - The DQ arrays of the average background and the target are combined
   using a bitwise-OR operation


If the target exposure is a simple ImageModel, the background image is
subtracted from it. If the target exposure is in the form of a CubeModel
(e.g. the result of a time-series exposure), the background image
is subtracted from all planes of the CubeModel.

The output results are always returned in a new data model, leaving the original
input model unchanged.

Upon successful completion of the step, the S_BKDSUB keyword will be set to
'COMPLETE' in the output product header.

Step Arguments
==============
There are no step-specific arguments.

Reference File
==============
The step does not use any reference files.
