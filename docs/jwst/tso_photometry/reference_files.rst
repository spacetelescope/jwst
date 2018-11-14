Reference File
==============
The tso_photometry step uses a TsoPhotModel reference file, reference type
TSOPHOT, that supplies values of radius (in pixels) for the target aperture
and the inner and outer radii for the background annulus.

.. include:: ../includes/standard_keywords.rst

.. include:: tsophot_selection.rst

Type Specific Keywords for TSOPHOT
----------------------------------
These keywords are required to be present in a TsoPhotModel reference file.

========  ====================
Keyword   Model Name
========  ====================
EXP_TYPE  meta.exposure.type
TSOVISIT  meta.visit.tsovisit
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
