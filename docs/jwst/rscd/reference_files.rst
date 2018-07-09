Reference File Types
====================

The RSCD correction step uses an RSCD reference file. 



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
It uses ``EXTNAME=RSCD`` and contains seventeen columns which are used in determining the correction
for the equations given after the table. 

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

In order to explain where these parameters are used in the correction we will go over the correction equations given
in the Description Section.

The general form of  the correction to be added to the input data is::

   corrected data = input data data + dn_accumulated * scale * exp(-T / tau)  (Equation 1)

where T is the time since the last group in the previous integration, tau is the exponential time constant found in the RSCD table
and  dn_accumulated is the DN level that was accumulated for the pixel from the previous integration.
In case where the last integration  does not saturate the :math:`scale` term in equation 1 is determined according to the equation:

       :math:`scale = b{1}* [Counts{2}^{b{2}} * [1/exp(Counts{2}/b{3}) -1 ]\; \; (Equation \; 2)`

The following two additional equations are used in Equation 2:

	  :math:`b{1} = ascale * (illum_{zpt} + illum_{slope}*N + illum2* N^2) \; \; (Equation \; 2.1)`
	  :math:`Counts{2} = Final \, DN \, in \, the \,  last \, group \, in \; the \, last \, integration 
	  \, - Crossover \, Point \; \; (Equation \; 2.2)`
	  
The parameters for equations 2, 2.1, and 2,2  are:
	  - :math:`b{2}` in equation 2 is table column POW from RSCD table
          - :math:`b{3}` in equation 2 is table column  PARAM3 from the RSCD table
	  - ascale  in equation 2.1 is in the RSCD table 
	  - :math:`illum_{zpt}`  in equation 2.1 is in the RSCD table 
	  - :math:`illum_{slope}`  in equation 2.1 is in the RSCD table
	  - :math:`illum2`  in equation 2.1 is in the RSCD table
	  - N  in equation 2.1 is the number of groups per integration
	  - Crossover Point in equation 2.2 is CROSSOPT in the RSCD table
	  
    
If the previous integration saturates, :math:`scale` is no longer calculated using equation 2 - 2.2, instead it is calculated 
using equations 3 and 3.1.

   :math:`scale_\text{sat} = slope * Counts{3} + sat_\text{mzp} \; \; (Equation \; 3)`
 
   :math:`slope = sat_{zp} + sat_{slope} * N + sat_2*N^2 + evenrow_{corrections} \; \; (Equation \; 3.1)`.

The parameters in equation 3 and 3.1 are:

    - :math:`Counts{3}`  in equation 3  is an estimate of the what the last group in the previous integration would 
    have been if saturation did not exist.
    - :math:`sat_\text{mzp}`  in equation 3 is in the RSCD table
    - :math:`scale_\text{sat}`  in equation 3 is SAT_SCALE in the RSCD table
    - :math:`sat_{zp}` in equation 3.1 is in the RSCD table
    - :math:`sat_{slope}` in equation 3.1 is  in the RSCD table
    - :math:`sat_2` in equation 3.1 is SAT2 in the RSCD table
    - :math:`evenrow_{corrections}` in equation 3.1 is SAT_ROWTERM in the RSCD table
    - N  is the number of groups per integration
 
   

