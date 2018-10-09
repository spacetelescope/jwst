Reference File
==============

There are two types of reference files used by the cube_build step. The first type holds the default
cube parameters used in setting up the output IFU Cube. The reftype for this reference file is *cubepar*
and there is a  reference file of this type for MIRI data and one for NIRSPEC data. These files contain tables 
for each band of the spatial and spectral size and the size of the region of interest to use to construct the
IFU cube.  If more than one band is used to build the IFU cube,
then the final spatial and spectral size will be the smallest one from the list of input bands. 
Currently cube_build can only produce IFU cubes with a linear spatial and spectral dimension. In the future we plan 
to allow a varying spectral step with wavelength.

The other type of reference file pertains only to MIRI data and contains the width of the PSF and LSF per
band. The reftype for this reference file is *resol*.
This information is used if the weight function incorporates the size of the psf and lsf, i.e.  --weighting = miripsf 

Cube Building Parameter Reference File Format
---------------------------------------------
The cube parameter reference files are FITS files with  BINTABLE extensions. The FITS primary data array is
assumed to be empty. The MIRI cube parameter file contains three  BINTABLE extensions, while the NIRSPEC  
file contains five BINTABLE extensions. In both files the first extension contains
the spatial and spectral cube sample size for each band. The second extension holds the  Modified Shepard weighting
values to use for each band. The third extension will be used in Build 7.2 and contains the  wavelengths and associated
region of interest size to use  if the IFU cubes are created from several bands and the final output is to
have an IFU cube of  varying spectral scale. In the case of MIRI  the twelve spectral
bands can be combined into a single IFU cube an all the information to create cubes of varying
wavelength sampling  is contained in this third BINTABLE extension.   However for NIRSPEC data there are 
three types of multi-band cubes: PRISM, MEDIUM and HIGH resolution.  The third, forth and fifth BINTABLE
extensions  in the NIRSPEC 
reference file contains the wavelength sampling and region of interest size  to use for 
PRISM, MEDIUM resolution, and HIGH resolution multi-band cubes, respectively.

.. include:: ../includes/standard_keywords.rst

.. include:: cubepar_selection.rst

.. include:: cubepar_format.rst

MIRI Resolution reference file
------------------------------
The MIRI resolution reference file is a FITS file with four BINTABLE extensions. The FITS primary data array is
assumed to be empty. The first  BINTABLE extension  contains the RESOLVING_POWER the information to use for 
each band. This table has 12 rows and 11 columns, one row of information for each band.  The parameters in the 11 columns
provide the polynomial coefficients to determine the resolving power for band that row corresponds to. 
The second BINTABLE extension, PSF_FWHM_ALPHA,
has a format of 1 row and 5 columns. The 5 columns hold the polynomial coefficients for determining the alpha PSF
size. 
The third BINTABLE extension, PSF_FWHM_BETA,
has a format of 1 row and 5 columns. The 5 columns hold the polynomial coefficients for determining the beta PSF
size. 

.. include:: resol_selection.rst

.. include:: resol_format.rst
