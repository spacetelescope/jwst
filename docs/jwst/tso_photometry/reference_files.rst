Reference File
==============
The tso_photometry step uses a TsoPhotModel reference file, reference type
TSOPHOT, that supplies values of radius (in pixels) for the target aperture
and the inner and outer radii for the background annulus.

CRDS Selection Criteria
-----------------------
TSOPHOT reference files are selected on the basis of INSTRUME, EXP_TYPE,
and TSOVISIT.  For MIRI exposures, EXP_TYPE should be MIR_IMAGE.  For
NIRCam exposures, EXP_TYPE should be NRC_TSIMAGE.  For both MIRI and NIRCam,
TSOVISIT should be True.

Required keywords
-----------------
These keywords are required to be present in a TsoPhotModel reference file.
The first column gives the FITS keyword names (although these reference
files are ASDF).  The second column gives the model name, which is needed
when creating and populating a new reference file.

========  ====================
Keyword   Model Name
--------  --------------------
AUTHOR    meta.author
DATAMODL  meta.model_type
DATE      meta.data
DESCRIP   meta.description
EXP_TYPE  meta.exposure.type
FILENAME  meta.filename
INSTRUME  meta.instrument.name
PEDIGREE  meta.pedigree
REFTYPE   meta.reftype
TELESCOP  meta.telescope
TSOVISIT  meta.visit.tsovisit
USEAFTER  meta.useafter
========  ====================

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
