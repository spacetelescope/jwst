Description
===========
At its basic level this step flat-fields an input science data set by dividing
by a flat-field reference image. In particular, the SCI array from the
flat-field reference file is divided into both the SCI and ERR arrays of the
science data set, and the flat-field DQ array is combined with the science DQ
array using a bit-wise OR operation.

MultiSlit data models, such as those resulting from the extracted regions of a
NIRSpec fixed slit observation, are handled as follows. First, if the
flat-field reference file supplied to the step is also in the form of a
MultiSlit model, it searches the reference file for slits with names that
match the slits in the science exposure (e.g. 'S1600A1' or 'S200B1'). When it
finds a match, it uses the flat-field data for that slit to correct the
particular slit data in the science exposure. If, on the other hand, the
flat-field consists of a single image model, the region corresponding to each
slit in the science data is extracted on-the-fly from the flat-field data and
applied to the corresponding slit in the science data.

Multiple-integration datasets (the _rateints.fits products from the ramp_fit
step) are handled by applying the flat-field to each integration.

For pixels whose DQ is NO_FLAT_FIELD in the reference file, the flat
value is reset to 1.0. Similarly, for pixels whose flat value is NaN, the flat
value is reset to 1.0 and DQ value in the output science data is set to
NO_FLAT_FIELD. In both cases, the effect is that no flat-field is applied.

If any part of the input data model gets flat-fielded (e.g. at least one
slit of a MultiSlit model), the status keyword S_FLAT will be set to
COMPLETE in the output science data.

Subarrays
---------
This step handles input science exposures that were taken in subarray modes in
a flexible way. If the reference data arrays are the same size as the science
data, they will be applied directly. If there is a mismatch, the routine will
extract a matching subarray from the reference file data arrays and apply them
to the science data. Hence full-frame reference files can be
used for both full-frame and subarray science exposures, or subarray-dependent
reference files can be provided if desired.
