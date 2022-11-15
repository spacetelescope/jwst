Description
===========

:Class: `jwst.flatfield.FlatFieldStep`
:Alias: flat_field

At its basic level this step flat-fields an input science dataset by dividing
by a flat-field reference image. In particular, the SCI array from the
flat-field reference file is divided into the SCI array of the
science dataset, the flat-field DQ array is combined with the science DQ
array using a bitwise OR operation, and variance and error arrays in the
science dataset are updated to include the flat-field uncertainty.
Details for particular modes are given in the sections below.

Upon completion of the step, the step status keyword "S_FLAT" gets set
to "COMPLETE" in the output science data.

Imaging and Non-NIRSpec Spectroscopic Data
------------------------------------------
Simple imaging data, usually in the form of an ImageModel, and some
spectroscopic modes, use a straight-forward approach that involves applying
a single flat-field reference file to the science image. The spectroscopic
modes included in this category are NIRCam WFSS and Time-Series Grism,
NIRISS WFSS and SOSS, and MIRI MRS and LRS. All of these modes are processed
as follows:

- If the science data have been taken using a subarray and the FLAT
  reference file is a full-frame image, extract the corresponding subarray
  region from the flat-field data.

- Find pixels that have a value of NaN or zero in the FLAT reference file
  SCI array and set their DQ values to "NO_FLAT_FIELD" and "DO_NOT_USE."

- Reset the values of pixels in the flat that have DQ="NO_FLAT_FIELD" to
  1.0, so that they have no effect when applied to the science data.

- Propagate the FLAT reference file DQ values into the science exposure
  DQ array using a bitwise OR operation.

- Apply the flat according to:

  .. math::
     SCI_{science} = SCI_{science} / SCI_{flat}

  .. math::
     VAR\_POISSON_{science} = VAR\_POISSON_{science} / SCI_{flat}^2

  .. math::
     VAR\_RNOISE_{science} = VAR\_RNOISE_{science} / SCI_{flat}^2

  .. math::
     VAR\_FLAT_{science} = ( SCI_{science}^{2} / SCI_{flat}^{2} ) * ERR_{flat}^{2}

  .. math::
     ERR_{science} = \sqrt{VAR\_POISSON + VAR\_RNOISE + VAR\_FLAT}

Multi-integration datasets ("_rateints.fits" products), which are common
for modes like NIRCam Time-Series Grism, NIRISS SOSS, and MIRI LRS Slitless,
are handled by applying the above equations to each integration.

For guider exposures, the flat is applied in the same manner as given
in the equations above, except for several differences.  First, the variances
due to Poisson noise and read noise are not calculated.  Second, the output
ERR array is the combined input ERR plus the flatfield ERR, summed in
quadrature.

NIRSpec Spectroscopic Data
--------------------------
Flat-fielding of NIRSpec spectrographic data differs from other modes
in that the flat-field array that will be applied to the science data
is not read directly from CRDS.  This is because the flat-field varies with
wavelength and the wavelength of light that falls on any given pixel
depends on the mode and which slits are open.  The flat-field array
that is divided into the SCI and ERR arrays is constructed on-the-fly
by extracting the relevant section from the reference files, and then --
for each pixel -- interpolating to the appropriate wavelength for that
pixel.  This interpolation requires knowledge of the dispersion direction,
which is read from keyword "DISPAXIS."  See the Reference File section for
further details.

For NIRSpec Fixed-Slit and MOS exposures, an on-the-fly flat-field is
constructed to match each of the slits/slitlets contained in the science
exposure. For NIRSpec IFU exposures, a single full-frame flat-field is
constructed, which is applied to the entire science image.

NIRSpec NRS_BRIGHTOBJ data are processed just like NIRSpec Fixed-Slit
data, except that NRS_BRIGHTOBJ data are stored in a CubeModel,
rather than a MultiSlitModel.  A 2-D flat-field image is constructed
on-the-fly as usual, but this image is then divided into each plane of
the 3-D science data arrays.

In all cases, there is a step option that allows for saving the
on-the-fly flatfield to a file, if desired.

NIRSpec Fixed-Slit Primary Slit
-------------------------------
The primary slit in a NIRSpec fixed-slit exposure receives special handling.
If the primary slit, as given by the "FXD_SLIT" keyword value, contains a
point source, as given by the "SRCTYPE" keyword, it is necessary to know the
flatfield conversion factors for both a point source and a uniform source
for use later in the :ref:`master background <master_background_step>` step
in Stage 3 processing. The point source version of the flatfield correction
is applied to the slit data, but that correction is not appropriate for the
background signal contained in the slit, and hence corrections must be
applied later in the :ref:`master background <master_background_step>` step.

So in this case the `flatfield` step will compute 2D arrays of conversion
factors that are appropriate for a uniform source and for a point source,
and store those correction factors in the "FLATFIELD_UN" and "FLATFIELD_PS"
extensions, respectively, of the output data product. The point source
correction array is also applied to the slit data.

Note that this special handling is only needed when the slit contains a
point source, because in that case corrections to the wavelength grid are
applied by the :ref:`wavecorr <wavecorr_step>` step to account for any
source mis-centering in the slit and the flatfield conversion factors are
wavelength-dependent. A uniform source does not require wavelength corrections
and hence the flatfield conversions will differ for point and uniform
sources. Any secondary slits that may be included in a fixed-slit exposure
do not have source centering information available, so the
:ref:`wavecorr <wavecorr_step>` step is not applied, and hence there's no
difference between the point source and uniform source flatfield
conversions for those slits.
