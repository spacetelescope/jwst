
Description
============
This step applies a fringe correction to the SCI data of an input data set
by dividing the SCI and ERR arrays by a fringe reference image. In particular,
the SCI array from the fringe reference file is divided into the SCI and ERR
arrays of the science data set. Only pixels that have valid values in the SCI
array of the reference file will be corrected.

This correction is applied only to MIRI MRS (IFU) mode exposures, which are
always single full-frame 2-D images.

The input to this step is always an ImageModel data model. The fringe reference
file that matches the input detector (MIRIFUSHORT or MIRIFULONG) and wavelength
band (SHORT, MEDIUM, or LONG, as specified by GRATNG14) is used.

Upon successful application of this correction, the status keyword S_FRINGE will
be set to COMPLETE.
