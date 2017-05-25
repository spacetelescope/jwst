Reference File Types
--------------------
The MIRI MRS straylight correction step uses a straylight mask. There
are three MIRI MRS SW masks, one for each of the three bands (SHORT,MEDIUM and LONG).

CRDS Selection Criteria
-----------------------
MIRI MRS reference files are selected on the basis of INSTRUME, DETECTOR, 
and BAND values from the input science data set.

MIRI MRS straylight Reference File Format
-----------------------------------------
The straylight mask  reference files are FITS files with one  
IMAGE extension. This image extension is a 2-D integer mask file of size 
1032 X 1024. The mask contains values of 1 for pixels that fall in 
the slice gaps and values of 0 for science pixels. The straylight 
algorithm only uses pixels that fall in the slice gaps to determine 

