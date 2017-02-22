Description
============

The photom step loads information into a data product that allows for the
conversion of count rates to absolute flux units. The flux conversion
information is read from the photometric reference file. The exact nature
of the information that's stored in the reference file and loaded into the
science data product depends on the instrument mode.

For most instrument modes the photom reference file contains a table of
exposure parameters that define various observing modes and the flux
conversion data for each of those modes. The table contains one row for each
allowed combination of exposure parameters,
such as detector, filter, pupil, and grating. The photom step searches the
table for the row that matches the parameters of the science exposure and
then copies the calibration information from that table row into the science
product.

For the MIRI MRS mode, the photom reference file contains arrays of sensitivity
factors and pixel sizes that are copied into the science product and also
applied to the SCI and ERR arrays of the science product.

For the table-based reference files, the calibration information in each row
includes a scalar conversion
constant, as well as optional arrays of wavelength and relative response
(as a function of wavelength). The scalar conversion constant in a selected
table row is copied into the keyword PHOTMJSR in the primary header of the
science product. The value of PHOTMJSR can then be used to convert data from
units of DN/sec to MJy/steradian. The step also computes, on the fly,
the equivalent conversion factor for converting the data to units of
microJy/square-arcsecond and stores this value in the header keyword PHOTUJA2.

If the photom step finds that the wavelength and relative response arrays are
populated in the selected table row, it copies those arrays to a table extension
called "RELSENS" in the science data product. For the MIRI MRS mode, the
sensitivity factors and pixel size arrays are multiplied together and copied to
an image extension called "RELSENS2D" in the science data product.

For multiple-integration datasets, the photom step can take either of these as 
input: a dataset containing the slope results for each integration in the 
exposure, or the dataset containing the single slope image that is the result 
of averaging over all integrations. 

Finally, if the science data are from an imaging mode, which is determined
from the value of the EXP_TYPE keyword, the data from a pixel area map
reference file will also be loaded into the science data product. The 2D
data array from the pixel area map will be copied into an image extension
called "AREA" in the science data product. For imaging mode exposures, the
values of the PIXAR_SR and PIXAR_A2 keywords in the photom reference table
will also be copied into keywords of the same name in the primary header of
the science data product.

Note that, except for MIRI MRS exposures, the pixel values of the science data
product are not actually changed by the photom step. All of the reference data
is simply attached to the science data product in some way.

Upon successful completion of this step, the status keyword S_PHOTOM will be
set to COMPLETE.
