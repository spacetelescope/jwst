Description
===========

:Class: `jwst.pathloss.PathlossStep`
:Alias: pathloss

Overview
--------
The ``pathloss`` step calculates and applies corrections that are
needed to account for various types of signal loss in spectroscopic data.
The motivation behind the correction is different for the MIRI, NIRSpec,
and NIRISS observing modes.
For MIRI LRS fixed slit, this correction simply accounts for light not
passing through the aperture.
For NIRSpec, this correction accounts aperture losses as well as losses
in the optical system due to light being scattered outside the grating.
For NIRISS SOSS data it corrects for the light that falls outside the
detector subarray in use.

Background
----------
The correction is applicable to MIRI LRS fixed-slit, NIRSpec IFU, MSA,
and FIXEDSLIT, and NIRISS SOSS data.
The description of how the NIRSpec reference files were created and how they are to be
applied to NIRSpec data is given in ESA-JWST-SCI-NRS-TN-2016-004 (P. Ferruit:
The correction of path losses for uniform and point sources).  The NIRISS algorithm
was provided by Kevin Volk.

Algorithm
---------

NIRSpec
+++++++
This step calculates a 1-D correction array as a function of wavelength by
interpolating in the pathloss reference file cube at the position of a point source target.
It creates 2 pairs of 1-D arrays, a wavelength array (calculated from the WCS applied to
the index of the plane in the wavelength direction) and a pathloss correction array
calculated by interpolating each plane of the pathloss cube at the position of
the source (which is taken from the datamodel).  Pairs of these arrays are computed
for both point source and uniform source data types.
For the uniform source pathloss calculation, there is no dependence on position
in the aperture/slit.

Once the 1-D correction arrays have been computed, both forms of the correction
(point and uniform) are interpolated, as a function of wavelength, into
the 2-D space of the slit or IFU data and attached to the output data model
(extensions "PATHLOSS_PS" and "PATHLOSS_UN") as a record of what was computed.
The form of the 2-D correction (point or uniform) that's appropriate for the
data is divided into the SCI and ERR arrays and propagated into the variance
arrays of the science data.

MIRI LRS
++++++++
The algorithm for MIRI LRS mode is largely the same as that for NIRSpec described
above, with the exception of the format in which the reference data are stored.
First, the position of the target on the detector is estimated from the target RA/Dec
given in the exposure header (TARG_RA, TARG_DEC keywords). This position is then
used to interpolate within the pathloss reference data to compute a 1-D pathloss
correction array. The 1-D pathloss correction is then interpolated into the 2-D
space of the data being corrected based on the wavelengths of each pixel in the
science data. The 2-D correction array is then applied to the science data and
stored (as a record of what was applied) in the output datamodel ("PATHLOSS_PS"
extension).

The MIRI LRS correction is only applicable to point source data. The step is
skipped if the SRCTYPE of the input data does not indicate a point source.

NIRISS SOSS
+++++++++++
The correction depends on column number in the science data and on the Pupil Wheel
position (keyword PWCPOS).  It is provided in the reference file as a FITS image of
3 dimensions (to be compatible with the NIRSpec reference file format).  The first
dimension is a dummy, while the second gives the dependence with row number, and the
third with Pupil Wheel position.  For the SUBSTEP96 subarray, the reference file
data has shape (1, 2040, 17).

The algorithm calculates the correction for each column by simply interpolating
along the Pupil Wheel position dimension of the reference file using linear
interpolation.  The 1-D vector of correction vs. column number is interpolated,
as a function of wavelength, into the 2-D space of the science image and divided
into the SCI and ERR arrays and propagated into the variance arrays.
The 2-D correction array is also attached to the datamodel (extension "PATHLOSS_PS")
as a record of what was applied.

Error Propagation
-----------------
As described above, the NIRSpec and NIRISS correction factors are divided into the
SCI and ERR arrays of the science data, and the square of the correction is divided
into the variance arrays (VAR_RNOISE, VAR_POISSON, VAR_FLAT) if they exist.
For MIRI LRS, the correction factors are multiplicative, hence they are multiplied
into the SCI and ERR arrays, and the square of the correction is multiplied into
the variance arrays.
