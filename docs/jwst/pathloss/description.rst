
Description
===========

Overview
--------

The pathloss correction step calculates the correction to apply to spectra
when the 1-d extraction is performed.  The motivation behind the correction
is different for NIRSPEC and NIRISS, while that for MIRI has not been
implemented yet.  For NIRSPEC, this correction accounts for losses
in the optical system due to light being scattered outside the grating, and
to light not passing through the aperture, while for NIRISS SOSS data it
corrects for the flux that falls outside the subarray.

Background
__________

The correction is applicable to NIRSPEC IFU, MSA and FIXEDSLIT exposure types,
to NIRISS SOSS data, and to MIRI LRS and MRS data, although the MIRI
correction has not been implemented yet.
The description of how the NIRSPEC reference files were created and how they are to be
applied to NIRSPEC data is given in ESA-JWST-SCI-NRS-TN-2016-004 (P. Ferruit:
The correction of path losses for uniform and point sources).  The NIRISS algorithm
was provided by Kevin Volk.

Algorithm
_________

NIRSPEC
-------

This step calculates the pathloss 1-d array as a function of wavelength by
interpolating in the pathloss cube at the position of the point source target.
It creates 2 pairs of 1-d arrays, a wavelength array (calculated from the WCS applied to
the index of the plane in the wavelength direction) and a pathloss array
calculated by interpolating each plane of the pathloss cube at the position of
the source (which is taken from datamodel).  There are pairs of these arrays for
both pointsource and uniformsource data types.

For the uniform source pathloss calculation, there is no dependence on position
in the aperture, so the array of pathlosses and calculated wavelengths are attached
to the datamodel.

Using 1-d arrays for the pathloss is different from what is suggested in the
Ferruit document, where it is recommended that 2-d arrays of pathloss correction are
attached to the data.  However, since the only variable in the 2-d array is the
wavelength, it was decided to simplify the process (and remove the possibility of
incorrect usage) by creating 1-d arrays of pathloss and wavelength, which are to
be applied at the time of 1-d extraction.

2-d arrays of the pathloss correction are also provided.

NIRISS
------

The correction depends on column number in the science data and on the Pupil Wheel
position (Keyword PWCPOS).  It is provided in the reference file as a FITS image of
3 dimensions (to be compatible with the NIRSPEC reference file format).  The first
dimension is a dummy, while the second gives the dependence with row number and the
third with Pupil Wheel position.  For the SUBSTEP96 subarray, the reference file
data has shape (1, 2040, 17).

The algorithm simply calculates the correction for each column by interpolating
along the Pupil Wheel position dimension of the reference file using linear
interpolation.  The 1-d vector of correction vs. column number is attached to the
science data in the PS extension, and can be obtained from the
ImageModel using the .pathloss_pointsource attribute.  This is a vector of length
2048 which gives the correction to be applied to each column of the science data,
in the sense that the correction should be divided into the data to correct it.
