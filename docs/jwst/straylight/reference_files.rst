Reference File Types
====================
The default algorithm  in the MIRI MRS stray-light correction step uses information contained
in the  meta data of the input image which maps each pixels to a slice or the  region between the
slices, also known as the slice gaps. This information was previously loaded from a reference file into the meta data by the assign_wcs
step. 
There is an option to use a more simplistic algorithm that uses  stray-light mask reference file.

CRDS Selection Criteria
-----------------------
If --method = "Nearest" option is used then the  MIRI MRS stray-light reference file is  selected on the basis of INSTRUME, DETECTOR, 
and BAND values of the input science data set.

MIRI MRS stray-light  Reference File Format
-------------------------------------------
The stray-light mask  reference files are FITS files with  and empty primary data
array and one IMAGE extension. This IMAGE extension is
a 2-D integer image  mask file of size 
1032 X 1024. The mask contains values of 1 for pixels that fall in 
the slice gaps and values of 0 for science pixels. The stray-light 
algorithm only uses pixels that fall in the slice gaps to determine 
the correction.

