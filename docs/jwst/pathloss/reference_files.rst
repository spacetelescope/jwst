Reference File
==============
The pathloss correction step uses a pathloss reference file.

.. include:: standard_keywords.rst

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
