Reference File
==============
The default MIRI MRS stray-light correction step uses a 'regions' reference file type. This information
is stored in the ASDF extension of the input_image. The information contained in the
'regions' has been read in during the assign_wcs step from the 
first extension in the MIRI Distortion file and stored in the ASDF extension of the input model.
The more simplistic algorithm (an older algorithm) uses the stray-light mask. There
are three MIRI MRS SW masks, one for each of the three bands (SHORT,MEDIUM and LONG).

CRDS Selection Criteria
-----------------------
The MIRI MRS stray-light reference files are selected on the basis of INSTRUME, DETECTOR, 
and BAND values of the input science data set.

MIRI MRS stray-light  Reference File Format
------------------------------------------
The regions file contains the 'slice number' of the science pixels and a value of zero for
the gap pixels. 
The stray-light mask  reference files are FITS files with  and empty primary data
array and one IMAGE extension. This IMAGE extension is
a 2-D integer image  mask file of size 
1032 X 1024. The mask contains values of 1 for pixels that fall in 
the slice gaps and values of 0 for science pixels. The stray-light 
algorithm only uses pixels that fall in the slice gaps to determine 


