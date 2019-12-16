Description
===========
The NIRSpec MSA imprint subtraction step removes patterns created in NIRSpec
MOS and IFU exposures by the MSA structure. This is accomplished by
subtracting a dedicated exposure taken with all MSA shutters closed and the
IFU entrance aperture blocked.

The step has two input parameters: the target exposure and the imprint
exposure. These arguments can be provided as either a file name
or a JWST data model.

The SCI data array of the imprint exposure is subtracted from the SCI array
of the target exposure. The DQ arrays of the two exposures are combined using
a bitwise logical OR operation. The ERR and variance arrays are not
currently used or modified.

Step Arguments
==============
The imprint subtraction step has no step-specific arguments.

Reference File
==============
The imprint subtraction step does not use any reference files.
