Reference File
==============
The pathloss correction step uses a pathloss reference file.

.. include:: ../includes/standard_keywords.rst

.. include:: pathloss_selection.rst
             
Only NIRSPEC IFU, FIXEDSLIT and MSA data, and NIRISS SOSS data perform a
pathloss correction.

Pathloss Reference File Format
------------------------------

The PATHLOSS reference files are FITS files with extensions for each
of the aperture types. The FITS primary data array is assumed to be empty.

The NIRSPEC IFU reference file just has four extensions, one pair for
point sources, and one pair for uniform sources.  In each pair, there are
either 3-d arrays for point sources, because the pathloss correction depends
on the position of the source in the aperture, or 1-d arrays for uniform
sources.  The pair of arrays are the pathloss correction itself as a function
of decenter in the aperture (pointsource only) and wavelength, and the variance
on this measurement (currently estimated).

The NIRSPEC FIXEDSLIT reference file has this FITS structure:

==== =======     ==========   ====== =============  ==========
No.    Name         Type      Cards   Dimensions    Format
==== =======     ==========   ====== =============  ==========
0    PRIMARY     PrimaryHDU      15   ()              
1    PS          ImageHDU        29   (21, 21, 21)  float64   
2    PSVAR       ImageHDU        29   (21, 21, 21)  float64   
3    UNI         ImageHDU        19   (21,)         float64   
4    UNIVAR      ImageHDU        19   (21,)         float64   
5    PS          ImageHDU        29   (21, 21, 21)  float64   
6    PSVAR       ImageHDU        29   (21, 21, 21)  float64   
7    UNI         ImageHDU        19   (21,)         float64   
8    UNIVAR      ImageHDU        19   (21,)         float64   
9    PS          ImageHDU        29   (21, 21, 21)  float64   
10   PSVAR       ImageHDU        29   (21, 21, 21)  float64   
11   UNI         ImageHDU        19   (21,)         float64   
12   UNIVAR      ImageHDU        19   (21,)         float64   
13   PS          ImageHDU        29   (21, 21, 21)  float64   
14   PSVAR       ImageHDU        29   (21, 21, 21)  float64   
15   UNI         ImageHDU        19   (21,)         float64   
16   UNIVAR      ImageHDU        19   (21,)         float64   
==== =======     ==========   ====== =============  ==========
 

HDU #1-4 are for the S200A1 aperture, while #5-8 are for S200A2,
#9-12 are for S200B1 and #13-16 are for S1600A1.  Currently there is
no information for the S400A1 aperture.

The NIRSPEC IFU reference file just has 4 extensions after the primary HDU,
as the behavious or each slice is considered identical.

The NIRSPEC MSASPEC reference file has 2 sets of 4 extensions, one for the 1x1
aperture size, and one for the 1x3 aperture size.  Currently there are no other
aperture sizes.

The NIRISS SOSS reference file has 1 extension in addition to the primary
header unit.  It contains a 3-dimensional array of float32 correction values.
The dimensions of the array are 1x2040x17.  The first dimension is a dummy to
force the array dimensionality to be the same as the NIRSPEC reference file
arrays.  The other 2 dimensions refer to the number of columns in the correction
(the same as the number of columns in the science data) and the range of
values for the Pupil Wheel position (PWCPOS).

The headers associated with the PS extensions should contain the WCS
information that describes what variables the correction depends on and
how they relate to the dimensions of the correction array.

For the NIRSPEC reference files (MSASPEC, FIXEDSLIT and IFU), the WCS keywords
should look like this:

======= ========== =========================================================================================
Keyword Value      Comment
======= ========== =========================================================================================
CRPIX1  1.0        Reference pixel in fastest dimension
CRVAL1  -0.5       Coordinate value at this reference pixel
CDELT1  0.05       Change in coordinate value for unit change in index
CTYPE1  'UNITLESS' Type of physical coordinate in this dimension
======= ========== =========================================================================================

This dimension expresses the decenter along the dispersion direction for a point source

======= ========== =========================================================================================
CRPIX2  1.0        Reference pixel in fastest dimension
CRVAL2  -0.5       Coordinate value at this reference pixel
CDELT2  0.05       Change in coordinate value for unit change in index
CTYPE2  'UNITLESS' Type of physical coordinate in this dimension
======= ========== =========================================================================================

This dimension expresses the decenter along the direction perpendicular to the dispersion for a point source

======= ========== =========================================================================================
CRPIX3  1.0        Reference pixel in fastest dimension
CRVAL3  6.0E-7     Coordinate value at this reference pixel
CDELT3  2.35E-7    Change in coordinate value for unit change in index
CTYPE3  'Meter'    Type of physical coordinate in this dimension (should be 'WAVELENGTH')
======= ========== =========================================================================================

This dimension expresses the change of correction with wavelength

The NIRISS SOSS reference file should also have WCS components, but their
interpretation is different from those in the NIRSPEC reference file:

======= ===================== =========================================================================================
Keyword Value                 Comment
======= ===================== =========================================================================================
CRPIX1  5.0                   Reference pixel in fastest dimension
CRVAL1  5.0                   Coordinate value at this reference pixel
CDELT1  1.0                   Change in coordinate value for unit change in index
CTYPE1  'PIXEL'               Type of physical coordinate in this dimension
======= ===================== =========================================================================================

This dimension expresses the decenter along the dispersion direction for a point source

======= ===================== =========================================================================================
CRPIX2  9.0                   Reference pixel in fastest dimension
CRVAL2  245.304               Coordinate value at this reference pixel
CDELT2  0.1                   Change in coordinate value for unit change in index
CTYPE2  'Pupil Wheel Setting' Type of physical coordinate in this dimension
======= ===================== =========================================================================================

This dimension expresses the decenter along the direction perpendicular to the dispersion for a point source

======= ===================== =========================================================================================
CRPIX3  1.0                   Reference pixel in fastest dimension
CRVAL3  1.0                   Coordinate value at this reference pixel
CDELT3  1.0                   Change in coordinate value for unit change in index
CTYPE3  'Dummy'               Type of physical coordinate in this dimension (should be 'WAVELENGTH')
======= ===================== =========================================================================================

This dimension expresses the change of correction with wavelength

