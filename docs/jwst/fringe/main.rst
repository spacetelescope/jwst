Description
============
The ``fringe`` step applies a fringe correction to MIRI MRS images.
In particular, the SCI array from a fringe reference file is divided into the
SCI and ERR arrays of the science data set. Only pixels that have valid (non-NaN)
values in the SCI array of the reference file will be corrected.
The DQ and variance arrays of the science exposure are not currently modified by
this step.

The input to this step is in the form of an ImageModel data model. The fringe reference
file that matches the input detector (MIRIFUSHORT or MIRIFULONG) and wavelength
band (SHORT, MEDIUM, or LONG, as specified by GRATNG14) is used.

Upon successful application of this correction, the status keyword "S_FRINGE" is
set to "COMPLETE".
