Description
===========

Overview
--------

The pathloss correction step calculates the correction to apply to spectra
when the 1-d extraction is performed.  This correction accounts for losses
in the optical system due to light being scattered outside the grating, and
to light not passing through the aperture.

Background
__________

The correction is applicable to NIRSPEC IFU, MSA and FIXEDSLIT exposure types,
to NIRISS SOSS data, and to MIRI LRS and MRS data, although the MIRI and NIRISS
corrections are not implemented in Build 7.
The description of how the reference files were created and how they are to be
applied to NIRSPEC data is given in ESA-JWST-SCI-NRS-TN-2016-004 (P. Ferruit:
The correction of path losses for uniform and point sources).

Algorithm
_________

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
