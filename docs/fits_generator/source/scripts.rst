Command-line scripts
====================

create_data directory
---------------------

create_data followed by a directory will process the proposal file
(generally a 5-digit string followed by '.prop') in that directory.
The proposal file contains the names of the FITS files to be processed
and the relationship between the exposures, allowing a unique
numbering scheme.

Each FITS file referred to in the exposure will be processed to make a
Level1b format JWST dataset with the pixel data flipped and/or rotated
to make it conform to the DMS coordinate system, in which all imaging
data has roughly the same orientation and parity on the sky.

The 5-digit string is used in the name of the Level 1b product, in that
file 12345.prop will make data of the form

jw12345aaabbb_cccdd_eeeee_DATATYPE_uncal.fits.

The numbers that fill in the other letter spaces come from the structure
of the proposal file, which is a sequence of nested levels.  As each
level is repeated, the number assigned to repesent that level increments
by 1.




