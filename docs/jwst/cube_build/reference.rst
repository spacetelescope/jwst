Reference File
==============

There are two types of reference files used by the cube_build step. The first type holds the default
cube parameters used in setting up the output IFU Cube.  There is a reference file of this type for
MIRI data and one for NIRSPEC data. These files contain a table for each band of the spatial and spectral 
size to use to construct the IFU Cube. For build 7.1, if more than one band is used to build the IFU cube,
then the final spatial and spectral size will be the smallest one from the list of input bands. 
The IFU cubes have a linear spatial and spectral dimension. In build 7.2 we plan to alway a varying spectral
size with wavelength.

The other type of reference file pertains only to MIRI data and contains the width of the psf and lsf per
band. This information is used if the weight function incorporates the size of the psf and lsf, i.e.  --weighting = miripsf 

CRDS Selection Criteria
-----------------------
The cube parameter reference file selection is based on Instrument. 


Cube Building  Reference File Format
------------------------------------------
Description coming 

