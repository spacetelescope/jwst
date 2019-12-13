Description
===========
At its basic level this step flat-fields an input science data set by dividing
by a flat-field reference image. In particular, the SCI array from the
flat-field reference file is divided into both the SCI and ERR arrays of the
science data set, and the flat-field DQ array is combined with the science DQ
array using a bitwise OR operation. Details for particular modes are given
in the sections below.

Imaging and Non-NIRSpec Spectroscopic Data
------------------------------------------
Simple imaging data, usually in the form of an ImageModel, is handled as
follows:

- Find pixels that have a value of NaN or zero in the FLAT reference file
  SCI array and set their DQ values to "NO_FLAT_FIELD."

- Reset the values of pixels in the flat that have DQ="NO_FLAT_FIELD" to
  1.0, so that they have no effect when applied to the science data.

- Apply the flat by dividing it into the science exposure SCI and ERR arrays.

- Propagate the FLAT reference file DQ values into the science exposure
  DQ array using a bitwise OR operation.

Spectroscopic data in the form of MultiSlit data models are handled as follows:

- If the flat-field reference file supplied to the step is also in the form of a
  MultiSlit model, search the reference file for slits with names that
  match the slits in the science exposure (e.g. 'S1600A1' or 'S200B1').

- When a match is found, use the flat-field data for that slit to correct the
  particular slit data in the science exposure, using the same procedures as
  outlined above for imaging data.

If, on the other hand, the flat-field consists of a single image model, the
region corresponding to each slit in the science data is extracted on-the-fly
from the flat-field data and applied to the corresponding slit in the science data.

Multi-integration datasets ("_rateints.fits" products) are handled by applying
the above flat-field procedures to each integration.

If any part of the input data model gets flat-fielded (e.g. at least one
slit of a MultiSlit model), the status keyword "S_FLAT" gets set to
"COMPLETE" in the output science data.

NIRSpec Spectroscopic Data
--------------------------
Flat-fielding of NIRSpec spectrographic data differs from other modes
in that the flat-field array that will be applied to the science data
is not read directly from CRDS.  This is because the flat-field varies with
wavelength, and the wavelength of light that falls on any given pixel
depends on mode and on which slit or slits are open.  The flat-field array
that is divided into the SCI and ERR arrays is constructed on-the-fly
by extracting the relevant section from the reference files, and then --
for each pixel -- interpolating to the appropriate wavelength for that
pixel.  This interpolation requires knowledge of the dispersion direction,
which is gotten from keyword "DISPAXIS."  See the Reference File section for
further details.  There is an option to save the on-the-fly flat field to
a file.

NIRSpec NRS_BRIGHTOBJ data are processed much like other NIRSpec
spectrographic data, except that NRS_BRIGHTOBJ data are in a CubeModel,
rather than a MultiSlitModel or ImageModel (used for IFU data).  A 2-D
flat field image is constructed on-the-fly as usual, but this image
is divided into each plane of the 3-D science data and error arrays,
resulting in an output CubeModel.

Subarrays
---------
This step handles input science exposures that were taken in subarray modes in
a flexible way. If the reference data arrays are the same size as the science
data, they will be applied directly. If there is a mismatch, the routine will
extract a matching subarray from the reference file data arrays and apply them
to the science data. Hence full-frame reference files can be
used for both full-frame and subarray science exposures, or subarray-dependent
reference files can be provided if desired.

Error Propagation
-----------------
The VAR_POISSON and VAR_RNOISE variance arrays are divided by the square
of the flat-field value for each pixel. A flat-field variance array,
VAR_FLAT, is created from the science exposure and flat-field reference
file data using the following formula:

.. math::
   VAR\_FLAT = ( SCI_{science}^{2} / SCI_{flat}^{2} ) * ERR_{flat}^{2}

The total ERR array is updated as the square root of the quadratic sum of
VAR_POISSON, VAR_RNOISE, and VAR_FLAT.
