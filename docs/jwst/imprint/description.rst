Description
===========

:Class: `jwst.imprint.ImprintStep`
:Alias: imprint

The NIRSpec MSA imprint subtraction step removes patterns created in NIRSpec
MOS and IFU exposures by the MSA structure. This is accomplished by
subtracting a dedicated exposure taken with all MSA shutters closed and the
IFU entrance aperture blocked.  These exposures are called "imprint" or "leakcal"
images.

The imprint subtraction step has two input parameters: the target exposure and a
list of one or more imprint/leakcal exposures. When called as a standalone step, these
arguments can be provided as either file names or JWST data models.  When the step is called
in the context of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline, imprint/leakcal file
names should be specified in the input `spec2` association file, labeled with the member
type of "imprint".

In the event that multiple imprint images are provided, the step uses the
metadata of the target and imprint exposures to find the imprint exposure
that matches the observation number (keyword "OBSERVTN"), dither pattern
position number (keyword "PATT_NUM"), and background target flag
(keyword "BKGDTARG") of the input exposure. The matching
imprint image is then subtracted from the target image. If no matching imprint
image is found, the step will be skipped, returning the input target image
unaltered.

When subtracting the imprint data model from the target data model,
the SCI data array of the imprint exposure is subtracted from the SCI array
of the target exposure, and the DQ arrays of the two exposures are combined using
a bitwise logical OR operation. The error and variance arrays are not
currently used or modified.

Step Arguments
==============
The imprint subtraction step has no step-specific arguments.

Reference File
==============
The imprint subtraction step does not use any reference files.
