Reference File
==============
The barshadow step does uses the barshadow reference file.


CRDS Selection Criteria
-----------------------
The Barshadow reference file is selected only for exposures with
EXP_TYPE=NRS_MSASPEC.  All other EXP_TYPEs should return N/A.

.. include:: ../includes/standard_keywords.rst

Reference File Format
---------------------

The barshadow reference file is a FITS file that contains four extensions:

======== ===== ============ =========
EXTNAME  NAXIS Dimensions   Data type
======== ===== ============ =========
DATA1X1  2     101x1001     float
VAR1X1   2     101x1001     float
DATA1X3  2     101x1001     float
VAR1X3   2     101x1001     float
======== ===== ============ =========

The slow direction has 1001 rows and gives the dependence of
the bar shadow correction on the Y location of a pixel from the center of the
shutter.  The fast direction has 101 rows and gives the dependence
of the bar shadow correction of wavelength.  The WCS keywords CRPIX1/2, CRVAL1/2
and CDELT1/2 tell how to convert the pixel location in the reference file into
a Y location and wavelength.  The initial version of the reference file has Y varying from -1.0
for row 1 to +1.0 at row 1001, and the wavelength varying from 0.6x10^-6m to 5.3x10^-6m.

The extension headers look like this:

======== = ==================== = ==========================
XTENSION = 'IMAGE   '           / Image extension           
BITPIX   =                  -64 / array data type           
NAXIS    =                    2 / number of array dimensions
NAXIS1   =                  101                             
NAXIS2   =                 1001                             
PCOUNT   =                    0 / number of parameters      
GCOUNT   =                    1 / number of groups          
EXTNAME  = 'DATA1x1 '           / extension name            
BSCALE   =                  1.0                             
BZERO    =                  0.0                             
BUNIT    = 'UNITLESS'                                       
CTYPE1   = 'METER   '                                       
CTYPE2   = 'UNITLESS'                                       
CDELT1   =              4.7E-08                             
CDELT2   =                0.002                             
CRPIX1   =                  1.0                             
CRPIX2   =                  1.0                             
CRVAL1   =                6E-07                             
CRVAL2   =                 -1.0                             
APERTURE = 'MOS1x1  '                                       
HEIGHT   =           0.00020161                             
======== = ==================== = ==========================

