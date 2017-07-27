Reference File
==============
There are three reference file types for the persistence step:
TRAPDENSITY, PERSAT, and TRAPPARS.

CRDS Selection Criteria
-----------------------
Persistence reference files are selected by INSTRUME and DETECTOR.

At the present time, there are no reference files for MIRI, and CRDS
will return "N/A" for the names of the files if the persistence step
is run on MIRI data, in which case the input will be returned unchanged
except that the primary header keyword S_PERSIS will will have been
set to 'SKIPPED'.

Reference File Formats
----------------------
The TRAPDENSITY reference file contains an IMAGE extension that gives
the density of traps at each pixel.

The PERSAT reference file contains an IMAGE extension that gives the
persistence saturation threshold (full well) at each pixel.

The TRAPPARS reference file contains a BINTABLE extension with four
float (double precision) columns:

* capture0: the coefficient of the exponential capture term
* capture1: minus the reciprocal of the capture e-folding time
* capture2: the "instantaneous" capture coefficient
* decay_param: minus the reciprocal of the decay e-folding time
