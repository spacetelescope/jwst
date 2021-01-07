Description
===========

Overview
--------
The ``pathloss`` step calculates and applies corrections that are
needed to account for various types of signal loss in spectroscopic data.
The motivation behind the correction
is different for NIRSpec and NIRISS, while that for MIRI has not been
implemented yet.  For NIRSpec, this correction accounts for losses
in the optical system due to light being scattered outside the grating, and
to light not passing through the aperture, while for NIRISS SOSS data it
corrects for the flux that falls outside the subarray.

Background
----------
The correction is applicable to NIRSpec IFU, MSA, and FIXEDSLIT exposure types,
to NIRISS SOSS data, and to MIRI LRS and MRS data, although the MIRI
correction has not been implemented yet.
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

NIRISS
++++++
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
As described above, the correction factors are divided into the SCI and ERR
arrays of the science data, and the square of the correction is divided into the
variance arrays (VAR_RNOISE, VAR_POISSON, VAR_FLAT) if they exist.
