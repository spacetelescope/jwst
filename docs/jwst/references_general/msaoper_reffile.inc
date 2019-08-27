.. _msaoper_reffile:

MSAOPER Reference File
----------------------

:REFTYPE: MSAOPER
:Data model: N/A

The MSAOPER reference file contains a list of failed MSA shutters and their
failure state (stuck open, stuck closed, etc.)

.. include:: ../references_general/msaoper_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for MSAOPER
++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in MSAOPER reference files,
because they are used as CRDS selectors
(see :ref:`msaoper_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
+++++++++++++++++++++
The MSAOPER reference files are json format.

The fields are:

:title: Short description of the reference file
:reftype: Should be "MSAOPER"
:pedigree: Should be one of "DUMMY", "GROUND" or "INFLIGHT"
:author: Creator of the file
:instrument: JWST Instrument; should be "NIRSPEC"
:exp_type: EXP_TYPEs this file should be used with; should be "NRS_IFU|NRS_MSASPEC"
:telescope: Should be "JWST"
:useafter: Exposure datetime after which this file is applicable
:descrip: Description of reference file
:msaoper:
    :Q: Quadrant; should be an integer 1-4
    :x: x location of shutter (integer, 1-indexed)
    :y: y location of shutter (integer, 1-indexed)
    :state: state of shutter; should be "closed" or "open"
    :TA state: TA state of shutter; should be "closed" or "open"
    :Internal state: Internal state of shutter; should be "closed", "normal" or "open"
    :Vignetted: Is the shutter vignetted?  Should be "yes" or "no"
:history: Description of the history relevant to this file; might point to documentation
