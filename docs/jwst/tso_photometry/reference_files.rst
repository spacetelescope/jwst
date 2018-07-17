Reference File
==============
The tso_photometry step uses a TSOPHOT reference file that supplies values
of radius (in pixels) for the target aperture and the inner and outer radii
for the background annulus.

CRDS Selection Criteria
-----------------------
TSOPHOT reference files are selected on the basis of INSTRUME, EXP_TYPE,
and TSOVISIT.  For MIRI exposures, EXP_TYPE should be MIR_IMAGE.  For
NIRCam exposures, EXP_TYPE should be NRC_TSIMAGE.  For both MIRI and NIRCam,
TSOVISIT should be True.

TSOPHOT Reference File Format
-----------------------------
TSOPHOT reference files are ASDF files.  An object called 'radii' in a
TSOPHOT file defines the radii that the step needs.  This object is a list
of one or more dictionaries.  Each such dictionary has four keys: 'pupil',
'radius', 'radius_inner', and 'radius_outer'.  The particular one of
these dictionaries to use is selected by comparing meta.instrument.pupil
with the value corresponding to 'pupil' in each dictionary.  If an exact
match is found, that dictionary will be used.  If no match is found, the
first dictionary with 'pupil': 'ANY' will be selected.  The radii will be
taken from the values of keys 'radius', 'radius_inner', and 'radius_outer'.
