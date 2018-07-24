Reference File
==============
The msaflagopen correction step uses a MSAOPER reference file.

CRDS Selection Criteria
:::::::::::::::::::::::

:NIRSPEC: USEAFTER

Msaoper reference files are selected on the basis of USEAFTER date only.
They are valid for NIRSpec only.

MSAOPER Reference File Format
:::::::::::::::::::::::::::::

The MSAOPER reference files are json files.

The fields are:

:title: Short description of the reference file
:reftype: Should be "MSAOPER"
:pedigree: Should be one of "DUMMY", "GROUND" or "INFLIGHT"
:author: Creator of the file
:instrument: JWST Instrument, should be "NIRSPEC"
:exp_type: EXP_TYPEs this file should be used with, should be "NRS_IFU|NRS_MSASPEC"
:telescope: Should be "JWST"
:useafter: Exposure datetime after which this file is applicable
:descrip: Description of reference file
:msaoper:
    :Q: Quadrant, should be an integer 1-4
    :x: x location of shutter (integer)
    :y: y location of shutter (integer)
    :state: state of shutter, should be "closed" or "open"
    :TA: state TA state of shutter, should be "closed" or "open"
    :Internal_state: Internal state of shutter, should be "closed", "normal" or "open"
    :Vignetted: Is the shutter vignetted?  Should be "yes" or "no"
:history: Description of the history relevant to this file, might point to documentation
