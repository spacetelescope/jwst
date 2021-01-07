.. _pathloss_reffile:

PATHLOSS Reference File
-----------------------

:REFTYPE: PATHLOSS
:Data model: `~jwst.datamodels.PathlossModel`

The PATHLOSS reference file contains correction factors as functions of
source position in the aperture and wavelength.

.. include:: ../references_general/pathloss_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for PATHLOSS
+++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in PATHLOSS reference files,
because they are used as CRDS selectors
(see :ref:`pathloss_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
+++++++++++++++++++++
The PATHLOSS reference files are FITS files with extensions for each
of the aperture types. The FITS primary HDU does not contain a data array.

NIRSpec IFU
~~~~~~~~~~~
The NIRSpec IFU PATHLOSS reference file just four extensions: one pair for
point sources and one pair for uniform sources.  In each pair, there are
either 3-D arrays for point sources, because the pathloss correction depends
on the position of the source in the aperture, or 1-D arrays for uniform
sources.  The pair of arrays are the pathloss correction itself as a function
of decenter in the aperture (pointsource only) and wavelength, and the variance
on this measurement (currently estimated). The data apply equally to all IFU
slices. The structure of the FITS file is as follows:

==== =======  ==========  =============  ==========
HDU  EXTNAME  XTENSION     Dimensions    Data type
==== =======  ==========  =============  ==========
1    PS       ImageHDU     21 x 21 x 21  float64   
2    PSVAR    ImageHDU     21 x 21 x 21  float64   
3    UNI      ImageHDU     21            float64   
4    UNIVAR   ImageHDU     21            float64   
==== =======  ==========  =============  ==========

NIRSpec Fixed Slit
~~~~~~~~~~~~~~~~~~
The NIRSpec Fixed Slit reference file has the following FITS structure:

==== ======= ====== ========== =============  ==========
HDU  EXTNAME EXTVER XTENSION    Dimensions    Data type
==== ======= ====== ========== =============  ==========
1    PS        1    ImageHDU    21 x 21 x 21  float64   
2    PSVAR     1    ImageHDU    21 x 21 x 21  float64   
3    UNI       1    ImageHDU    21            float64   
4    UNIVAR    1    ImageHDU    21            float64   
5    PS        2    ImageHDU    21 x 21 x 21  float64   
6    PSVAR     2    ImageHDU    21 x 21 x 21  float64   
7    UNI       2    ImageHDU    21            float64   
8    UNIVAR    2    ImageHDU    21            float64   
9    PS        3    ImageHDU    21 x 21 x 21  float64   
10   PSVAR     3    ImageHDU    21 x 21 x 21  float64   
11   UNI       3    ImageHDU    21            float64   
12   UNIVAR    3    ImageHDU    21            float64   
13   PS        4    ImageHDU    21 x 21 x 21  float64   
14   PSVAR     4    ImageHDU    21 x 21 x 21  float64   
15   UNI       4    ImageHDU    21            float64   
16   UNIVAR    4    ImageHDU    21            float64   
==== ======= ====== ========== =============  ==========
 
HDU's 1--4 are for the S200A1 aperture, 5--8 are for S200A2,
9--12 are for S200B1, and 13--16 are for S1600A1.  Currently there is
no reference data for the S400A1 aperture.

NIRSpec MSASPEC
~~~~~~~~~~~~~~~
The NIRSpec MSASPEC reference file has 2 sets of 4 extensions: one set for the 1x1
aperture (slitlet) size and one set for the 1x3 aperture (slitlet) size.
Currently there is not any reference data for other aperture sizes.
The FITS file has the following structure:

==== ======= ====== ========== =============  ==========
HDU  EXTNAME EXTVER XTENSION    Dimensions    Data type
==== ======= ====== ========== =============  ==========
1    PS        1    ImageHDU    21 x 63 x 21  float64   
2    PSVAR     1    ImageHDU    21 x 63 x 21  float64   
3    UNI       1    ImageHDU    21            float64   
4    UNIVAR    1    ImageHDU    21            float64   
5    PS        2    ImageHDU    21 x 63 x 21  float64   
6    PSVAR     2    ImageHDU    21 x 63 x 21  float64   
7    UNI       2    ImageHDU    21            float64   
8    UNIVAR    2    ImageHDU    21            float64   
==== ======= ====== ========== =============  ==========

NIRISS SOSS
~~~~~~~~~~~
The NIRISS SOSS reference file has just 1 extension HDU.
The structure of the FITS file is as follows:

==== ======= ========== =============  ==========
HDU  EXTNAME XTENSION    Dimensions    Data type
==== ======= ========== =============  ==========
1    PS      ImageHDU   17 x 2040 x 1  float32   
==== ======= ========== =============  ==========

The PS extension contains a 3-D array of correction values.
The third dimension (length = 1) is a dummy to
force the array dimensionality to be the same as the NIRSpec reference file
arrays.  The other 2 dimensions refer to the number of columns in the correction
(the same as the number of columns in the science data) and the range of
values for the Pupil Wheel position (PWCPOS).

WCS Header Keywords
~~~~~~~~~~~~~~~~~~~
The headers of the PS extensions in all of the above reference
files should contain the WCS information that describes what variables the
correction depends on and how they relate to the dimensions of the correction array.

For the NIRSpec reference files (IFU, Fixed Slit, and MSASPEC), the WCS keywords
should have the values shown in the tables below.
Dimension 1 expresses the decenter along the dispersion direction for a point source:

======= ========== ==============================================================================
Keyword Value      Comment
======= ========== ==============================================================================
CRPIX1  1.0        Reference pixel in fastest dimension
CRVAL1  -0.5       Coordinate value at this reference pixel
CDELT1  0.05       Change in coordinate value for unit change in index
CTYPE1  'UNITLESS' Type of physical coordinate in this dimension
======= ========== ==============================================================================

Dimension 2 expresses the decenter along the direction perpendicular to the dispersion
for a point source:

======= ========== ==============================================================================
CRPIX2  1.0        Reference pixel in fastest dimension
CRVAL2  -0.5       Coordinate value at this reference pixel
CDELT2  0.05       Change in coordinate value for unit change in index
CTYPE2  'UNITLESS' Type of physical coordinate in this dimension
======= ========== ==============================================================================

Dimension 3 expresses the change of correction as a function of wavelength:

======= ========== ==============================================================================
CRPIX3  1.0        Reference pixel in fastest dimension
CRVAL3  6.0E-7     Coordinate value at this reference pixel
CDELT3  2.35E-7    Change in coordinate value for unit change in index
CTYPE3  'METER'    Type of physical coordinate in this dimension
======= ========== ==============================================================================

The NIRISS SOSS reference file should also have WCS components, but their
interpretation is different from those in the NIRSpec reference file.
Dimension 1 expresses the column number in the science data:

======= ===================== ===================================================================
Keyword Value                 Comment
======= ===================== ===================================================================
CRPIX1  5.0                   Reference pixel in fastest dimension
CRVAL1  5.0                   Coordinate value at this reference pixel
CDELT1  1.0                   Change in coordinate value for unit change in index
CTYPE1  'PIXEL'               Type of physical coordinate in this dimension
======= ===================== ===================================================================

Dimension 2 expresses the value of the PWCPOS keyword:

======= ===================== ===================================================================
CRPIX2  9.0                   Reference pixel in fastest dimension
CRVAL2  245.304               Coordinate value at this reference pixel
CDELT2  0.1                   Change in coordinate value for unit change in index
CTYPE2  'Pupil Wheel Setting' Type of physical coordinate in this dimension
======= ===================== ===================================================================

Dimension 3 is a dummy axis for the NIRISS SOSS reference file:

======= ================= ======================================================================
CRPIX3  1.0               Reference pixel in fastest dimension
CRVAL3  1.0               Coordinate value at this reference pixel
CDELT3  1.0               Change in coordinate value for unit change in index
CTYPE3  'Dummy'           Type of physical coordinate in this dimension
======= ================= ======================================================================

