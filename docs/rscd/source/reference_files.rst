Reference File
==============
The RSCD correction step uses an RSCD reference file.  The correction
to be added to the input data has the form::

    corrected = input + dn_accumulated * scale * exp(-T / tau)

CRDS Selection Criteria
-----------------------
RSCD reference files are selected on the basis of INSTRUME and DETECTOR
values for the input science data set.  The row of the table in the
reference file is selected on the basis of READPATT and SUBARRAY.

RSCD Reference File Format
---------------------------
The RSCD reference files are FITS files with a BINTABLE extension. The FITS
primary data array is assumed to be empty.

The BINTABLE extension contains the row-selection criterea (subarray and
readpatt) and the parameters for a double-exponential correction function.
It uses ``EXTNAME=RSCD`` and contains six columns:

* subarray: string, FULL or a subarray name
* readpatt: string, SLOW or FAST
* tau1: float, e-folding time scale for the first exponential (unit is frames)
* scale1: float, scale factor for the first exponential
* tau2: float, e-folding time scale for the second exponential (frames)
* scale2: float, scale factor for the second exponential
