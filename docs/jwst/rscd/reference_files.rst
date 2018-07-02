Reference File Types
====================

The RSCD correction step uses an RSCD reference file. 


Description
-----------

This correction is only applied to integration > 1. 
The correction to be added to the input data has the form::

    corrected data = input data data + dn_accumulated * scale * exp(-T / tau)

where, dn_accumulated is the DN level that was accumulated for the pixel from the previous integration. 
In case where the dn_accumulated does not saturate the scale factor is determined as follows:

       :math:`scale = b{1}* [Counts{2}^{b{2}} * [1/exp(Counts{2}/b{3}) -1]`
    
where :math:`b{1} = ascale * (illum_{zpt} + illum_{slope}*N + illum2* N^2)` (N is the number of groups per integration)
and :math:`Counts{2}` = Final DN in the last frame in the last integration - Crossover Point (found in the reference file).

If the previous integration saturates, :math:`scale` is modified in the following manner:

   :math:`scale_\text{sat} = slope * Counts{3} + sat_\text{mzp}`, 
   
where :math:`Counts{3}` = extrapolated counts past saturation * sat_scale and :math:`slope = sat_{zp} + sat_{slope} * N + sat_2*N^2 + evenrow_{corrections}`.


All fourteen  parameters :math:`tau`, :math:`b{1}`, :math:`b{2}`, :math:`b{3}`, :math:`illum_{zpt}`,
:math:`illum_{slope}`, :math:`illum2`, :math:`Crossover Point`, :math:`sat_{zp}`, :math:`sat_{slope}`, :math:`sat_2`,
:math:`sat_{scale}`, :math:`sat_\text{mzp}`, and :math:`evenrow_{corrections}` are found in the RSCD reference files.

CRDS Selection Criteria
-----------------------
RSCD reference files are selected on the basis of INSTRUME and DETECTOR
values for the input science data set.  The reference file for each detector is a table of values based on
READPATT (FAST, SLOW) , SUBARRAY (FULL or one the various subarray types) , and ROWS type (even or odd row).
The fourteen correction values are read in separately for even and odd rows for the readout pattern and  
if it is for the full array or one of the imager subarrays. 

RSCD Reference File Format
---------------------------
The RSCD reference files are FITS files with a BINTABLE extension. The FITS
primary data array is assumed to be empty.

The BINTABLE extension contains the row-selection criteria (SUBARRAY, READPATT, and ROW type)  
and the parameters for a double-exponential correction function.
It uses ``EXTNAME=RSCD`` and contains seventeen columns:

* SUBARRAY: string, FULL or a subarray name
* READPATT: string, SLOW or FAST
* ROWS: string, EVEN or ODD
* TAU: float, e-folding time scale for the first exponential (unit is frames)
* ASCALE: float,  b1 in equation 
* POW: float, b2 in equation
* ILLUM_ZP: float
* ILLUM_SLOPE: float
* ILLUM2: float
* PARAM3: :math:`b{3}` in equation
* CROSSOPT: float, crossover point
* SAT_ZP: float
* SAT_SLOPE: float
* SAT2: float
* SAT_MZP: float
* SAT_ROWTERM: float
* SAT_SCALE: float

