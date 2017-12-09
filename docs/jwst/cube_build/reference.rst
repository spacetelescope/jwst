Reference File
==============

There are two types of reference files used by the cube_build step. The first type holds the default
cube parameters used in setting up the output IFU Cube. The reftype for this reference file is *cubepars*
and there is a  reference file of this type for MIRI data and one for NIRSPEC data. These files contain a table 
for each band of the spatial and spectral 
size to use to construct the IFU Cube. For build 7.1, if more than one band is used to build the IFU cube,
then the final spatial and spectral size will be the smallest one from the list of input bands. 
The IFU cubes have a linear spatial and spectral dimension. In build 7.2 we plan to alway a varying spectral
size with wavelength.

The other type of reference file pertains only to MIRI data and contains the width of the PSF and LSF per
band. The reftype for this reference file is *resol*.
This information is used if the weight function incorporates the size of the psf and lsf, i.e.  --weighting = miripsf 


CRDS Selection Criteria
-----------------------
The cube parameter reference file selection is based on Instrument. CRDS selection criteria for the MIRI resolution 
reference file is  also based on Instrument (a N/Q is returned for NIRSPEC data).


Cube Building Parameter Reference File Format
---------------------------------------------
The cube parameter reference files are FITS file with a BINTABLE extension. The FITS primary data array is
assumed to be empty. The BINTABLE extension contains the information on the default sample sizes of spatial
and spectral dimension of the output spectral cubes, as well as the size region of interest to use 
around each spaxel in chosing the detetector pixels to combine. 
The first two colunms in  reference files defines the which band the row describes. For MIRI column 1 holds
the channel and column two holds the sub-channel, while for NIRSPEC column 1 holds the grating and column 2
contains the filter. For each band defined by columns 1 and 2, columns 3-6 contain the  
spatial size of the output cube, the spectral size of the output cube, spatial region of interest size,
and wavelength region of interest. 


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

