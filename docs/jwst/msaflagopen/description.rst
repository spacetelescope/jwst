Description
===========

:Class: `jwst.msaflagopen.MSAFlagOpenStep`
:Alias: msa_flagging

Overview
--------
The ``msaflagopen`` step flags pixels in NIRSpec exposures that are affected by
MSA shutters that are stuck in the open position.

Background
----------
The correction is applicable to NIRSpec IFU and MSA exposure types.

Algorithm
---------
The set of shutters whose state is not commandable (i.e. they are permanently stuck
in 'open' or 'closed' positions) is recorded in the MSAOPER reference file.
The reference file is searched for all shutters with any of the quantities
'Internal state', 'TA state' or 'state' set to 'open'.

The step loops over the list of stuck open shutters.  For each shutter, the bounding box
that encloses the projection of the shutter onto the detector array is calculated,
and for each pixel in the bounding box, the WCS is calculated.  If the pixel is inside
the region affected by light through the shutter, the WCS will have valid values,
whereas if the pixel is outside, the WCS values will be NaN.  The indices of each non-NaN
pixel in the WCS are used to alter the corresponding pixels in the DQ array by OR'ing
their DQ value with that for "MSA_FAILED_OPEN."
