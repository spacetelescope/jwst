.. _tsophot_reffile:

TSOPHOT Reference File
----------------------

:REFTYPE: TSOPHOT
:Data model: `~jwst.datamodels.TsoPhotModel`

The TSOPHOT reference file contains photometry aperture radii for
TSO imaging observation modes.

.. include:: ../references_general/tsophot_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for TSOPHOT
++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in TSOPHOT reference files,
because they are used as CRDS selectors
(see :ref:`tsophot_selectors`):

=========  =========================
Keyword    Data Model Name
=========  =========================
EXP_TYPE   model.meta.exposure.type
TSOVISIT   model.meta.visit.tsovisit
=========  =========================

Reference File Format
+++++++++++++++++++++
TSOPHOT reference files are ASDF format.
An object called 'radii' in a
TSOPHOT file defines the radii that the step needs.  This object is a list
of one or more dictionaries.  Each such dictionary has four keys: 'pupil',
'radius', 'radius_inner', and 'radius_outer'.  The particular one of
these dictionaries to use is selected by comparing meta.instrument.pupil
with the value corresponding to 'pupil' in each dictionary.  If an exact
match is found, that dictionary will be used.  If no match is found, the
first dictionary with 'pupil': 'ANY' will be selected.  The radii will be
taken from the values of keys 'radius', 'radius_inner', and 'radius_outer'.
