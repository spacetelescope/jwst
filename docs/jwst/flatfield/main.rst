Description
===========
At its basic level this step flat-fields an input science data set by dividing
by a flat-field reference image. In particular, the SCI array from the
flat-field reference file is divided into both the SCI and ERR arrays of the
science data set, and the flat-field DQ array is combined with the science DQ
array using a bit-wise OR operation.

Non-NIRSpec Data
----------------
MultiSlit data models are handled as follows. First, if the
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

NIRSpec imaging data are corrected the same as non-NIRSpec data,
i.e. they will just be divided by a flat-field reference image.

For pixels whose DQ is NO_FLAT_FIELD in the reference file, the flat
value is reset to 1.0. Similarly, for pixels whose flat value is NaN, the flat
value is reset to 1.0 and DQ value in the output science data is set to
NO_FLAT_FIELD. In both cases, the effect is that no flat-field is applied.

If any part of the input data model gets flat-fielded (e.g. at least one
slit of a MultiSlit model), the status keyword S_FLAT will be set to
COMPLETE in the output science data.

NIRSpec Spectroscopic Data
--------------------------
Flat-fielding of NIRSpec spectrographic data differs from other modes
in that the flat field array that will be
divided into the SCI and ERR arrays of the input science data set is not
read directly from CRDS.  This is because the flat field varies with
wavelength, and the wavelength of light that falls on any given pixel
depends on mode and on which slit or slits are open.  The flat-field array
that is divided into the SCI and ERR arrays is constructed on-the-fly
by extracting the relevant section from the reference files, and then --
for each pixel -- interpolating to the appropriate wavelength for that
pixel.  This interpolation requires knowledge of the dispersion direction,
which is gotten from keyword DISPAXIS.  See the Reference File section for
further details.  There is an option to save the on-the-fly flat field to
a file.

NIRSpec NRS_BRIGHTOBJ data are processed much like other NIRSpec
spectrographic data, except that NRS_BRIGHTOBJ data are in a CubeModel,
rather than a MultiSlitModel or ImageModel (used for IFU data).  A 2-D
flat field image will be constructed on-the-fly as usual, but this image
will be divided into each plane of the 3-D science data and error array,
resulting in an output CubeModel.

When this step is called with NIRSpec imaging data as input, the data will be
flat-fielded as described in the section for non-NIRSpec data.

Subarrays
---------
This step handles input science exposures that were taken in subarray modes in
a flexible way. If the reference data arrays are the same size as the science
data, they will be applied directly. If there is a mismatch, the routine will
extract a matching subarray from the reference file data arrays and apply them
to the science data. Hence full-frame reference files can be
used for both full-frame and subarray science exposures, or subarray-dependent
reference files can be provided if desired.
