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
For fixed slit data, if the ``wavecorr`` step has been run to provide wavelength
corrections to point sources, the corrected wavelengths will be used to
calculate the point source pathloss, whereas the uncorrected wavelengths (appropriate 
for uniform sources) will be used to calculate the uniform source pathlosses.
The form of the 2-D correction (point or uniform) that's appropriate for the
data is divided into the SCI and ERR arrays and propagated into the variance
arrays of the science data.

The MSA reference file contains 2 entries: one for a 1x1 slit and one for a 1x3 slit.
Each entry contains the pathloss correction for point source and uniform sources.
The former depends on the position of the target in the fiducial shutter and
wavelength, whereas the latter depends on wavelength only.  The point source 
entry consists of a 3-d array, where 2 of the dimensions map to the location
of the source (ranging from -0.5 to 0.5 in both X and Y), while the third dimension
carries the wavelength dependence.  The 1x3 shutter is 3 times as large in Y as in X.

The entry to use for a point source target is determined by looking at the shutter_state
attribute of the slit used.  This is a string with a length equal to the number
of shutters that make up the slit, with 1 denoting an open shutter, 0 a closed
shutter and x the fiducial (target) shutter.  The reference entry is determined
by how many shutters next to the fiducial shutter are open:

If both adjacent shutters are closed, the 1x1 entry is used.  A matching
shutter_state might be 'x' or '10x01'

If both adjacent shutters are open, the center region of the 1x3 entry is used.
This would be the case for a slit with shutter state '1x1' or '1011x1'.

If one adjacent shutter is open and one closed, the 1x3 entry is used.  If the
shutter below the fiducial is open and the shutter above closed, then the upper
region of the 1x3 pathloss array is used.  This is implemented by adding 1 to the
Y coordinate of the target position (bringing it into the range +0.5 to +1.5),
moving it to the upper third of the pathloss array.  A matching shutter state
might be '1x' or '11x011'

Similarly, if the shutter below the fiducial is closed and that above is open, the
lower third of the pathloss array is used by subtracting 1 from the Y coordinate of
the target position (bringing it into the range -1.5 to -0.5).  A matching shutter
state could be 'x111' or '110x1'.

Once the X and Y coordinates of the source are mapped into a pixel location in the
spatial dimensions of the pathloss array using the WCS of the transformation of position
to pixel location, the wavelength dependence is determined
by interpolating at that (fractional) pixel position in each wavelength plane,
resulting in a pair of 1-d arrays of pathloss correction and wavelength.  These arrays
are used to interpolate the correction for each pixel of the 2-d extracted science
array, since each pixel has a different wavelength, and the correction is applied
to the science pixel array.

For uniform sources, there is no dependence of the pathloss correction on position,
so the correction arrays are just 1-d arrays of correction and wavelength.  The
correction depends only on the number of shutters in the slit:

If there is 1 shutter, the 1x1 entry is used

If there are 3 or more shutters, the 1x3 entry is used

If there are 2 shutters, the correction used is the average of the 1x1
and 1x3 entries.

Like for the point source case, the 1-d arrays of pathloss correction and wavelength
are used to interpolate the correction for each pixel in the science data, using the
wavelength of each pixel to interpolate into the pathloss correction array.

Upon successful completion of the step, the status keyword "S_PTHLOS"
in the primary header is set to "COMPLETE".  For each SCI extension, the "PTHLOSS"
keyword is set to "POINT" if the point source correction type was applied to the
slit data. It is set to "UNIFORM" if the uniform correction type was applied.
The keyword is not set if no correction was applied.

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

If for any reason the source is determined to be outside of the slit, the
correction defaults to the center of the slit.

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
