Description
===========

:Class: `jwst.fringe.fringe_step.FringeStep`
:Alias: fringe

The ``fringe`` step applies a fringe correction to MIRI MRS images.
In particular, the SCI array from a :ref:`fringe_reffile` is divided into the
SCI and ERR arrays of the science data set. The variance arrays are divided by the square of the reference data.
Only pixels that have valid (non-NaN)
values in the SCI array of the reference file will be corrected.
The DQ of the science exposure is not currently modified by
this step.

The input to this step is in the form of an `~stdatamodels.jwst.datamodels.ImageModel` data model. The fringe reference
file that matches the input detector (MIRIFUSHORT or MIRIFULONG) and wavelength
band (SHORT, MEDIUM, or LONG, as specified by GRATNG14) is used.

Upon successful application of this correction, the status keyword "S_FRINGE" is
set to "COMPLETE".

Step Arguments
--------------
The ``fringe`` step has no step-specific arguments.
