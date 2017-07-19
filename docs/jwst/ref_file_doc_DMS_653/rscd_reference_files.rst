Reference File Types
--------------------
The RSCD correction step uses an RSCD reference file. This correction only applied to 
integrations > 1.  The correction
to be added to the input data has the form::

    corrected = input + dn_accumulated * scale * exp(-T / tau)

where dn_accumulated is the DN level that was accumulated for the pixel from the
previous integration. 


CRDS Selection Criteria
-----------------------
RSCD reference files are selected on the basis of INSTRUME and DETECTOR
values for the input science data set.  The reference file for each detector 
is a table of values based on READPATT (FAST, SLOW), SUBARRAY (FULL or one the 
various subarray types), and ROWS type (even or odd row).  The correction values, 
tau and scale, are read in separately for even and odd rows, based on the 
readout pattern and if it is for the full array or one of the imager subarrays. 
The table actually contains the parameters for a double-exponential function, 
but currently only the single exponential values are used. 


RSCD Reference File Format
---------------------------
The RSCD reference files are FITS files with a BINTABLE extension. The FITS
primary data array is assumed to be empty.

The BINTABLE extension contains the row-selection criteria (SUBARRAY,
READPATT, and ROW type) and the parameters for a double-exponential 
correction function.  It uses ``EXTNAME=RSCD`` and contains seven columns:

* SUBARRAY: string, FULL or a subarray name
* READPATT: string, SLOW or FAST
* ROWS: string, EVEN or ODD
* TAU1: float, e-folding time scale for the first exponential (unit is frames)
* SCALE1: float, scale factor for the first exponential
* TAU2: float, e-folding time scale for the second exponential (frames)
* SCALE2: float, scale factor for the second exponential
