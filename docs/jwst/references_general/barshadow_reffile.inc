.. _barshadow_reffile:

BARSHADOW Reference File
------------------------

:REFTYPE: BARSHADOW
:Data model: `~jwst.datamodels.BarshadowModel`

.. include:: ../references_general/barshadow_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for BARSHADOW
++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in BARSHADOW reference files,
because they are used as CRDS selectors
(see :ref:`barshadow_selectors`):

=========  ========================
Keyword    Data Model Name
=========  ========================
EXP_TYPE   model.meta.exposure.type
=========  ========================

Reference File Format
+++++++++++++++++++++
BARSHADOW reference files are FITS format, with 4 IMAGE extensions.
The FITS primary data array is assumed to be empty.
The format and content of the file is as follows:

======== ======== ===== ============ =========
EXTNAME  XTENSION NAXIS Dimensions   Data type
======== ======== ===== ============ =========
DATA1X1  IMAGE      2     101x1001     float
VAR1X1   IMAGE      2     101x1001     float
DATA1X3  IMAGE      2     101x1001     float
VAR1X3   IMAGE      2     101x1001     float
======== ======== ===== ============ =========

The slow direction has 1001 rows and gives the dependence of the bar shadow
correction on the Y location of a pixel from the center of the shutter.  The
fast direction has 101 rows and gives the dependence of the bar shadow
correction of wavelength.  The WCS keywords CRPIX1/2, CRVAL1/2 and CDELT1/2
tell how to convert the pixel location in the reference file into a Y location
and wavelength.  The initial version of the reference file has Y varying from
-1.0 for row 1 to +1.0 at row 1001, and the wavelength varying from
0.6x10\ :sup:`-6` m to 5.3x10\ :sup:`-6` m.

An example extension header is as follows:

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

