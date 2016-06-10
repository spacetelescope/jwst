Description
============

The photom step loads information into a data product that allows for the
conversion of count rates to absolute flux units. The flux conversion
information is read from the photometric reference table. The exact nature
of the information that's stored in the reference table and loaded into the
science data product depends on the observing mode.

The photom reference table contains one row of information for each
combination of exposure parameters that defines a particular observing mode,
such as detector, filter, pupil, and grating. The photom step searches the
reference table for the row that matches the parameters of the exposure and
then copies the information from that table row into the science product.

The flux conversion information in each table row includes a scalar conversion
constant, as well as optional arrays of wavelength and relative response
(as a function of wavelength). The scalar conversion constant in a selected
table row is copied into the keyword PHOTMJSR in the primary header of the
science product. The value of PHOTMJSR can then be used to convert data from
units of counts/sec to MJy/steradian. The step also computes, on the fly,
the equivalent conversion factor for converting the data to units of
microJy/square-arcsecond and stores this value in the header keyword PHOTUJA2.

If the photom step finds that the wavelength and relative response arrays are
populated in the selected row, it copies those arrays to a table extension
called "RELSENS" in the science data product.

For multiple-integration datasets, the photom step can take either of these as 
input: a dataset containing the slope results for each integration in the 
exposure, or the dataset containing the single slope image that is the result 
of averaging over all integrations. 

Finally, if the science data are from an imaging mode, which is determined
from the value of the EXP_TYPE keyword, the data from a pixel area map
reference file will also be loaded into the science data product. The 2-D
data array from the pixel area map will be copied into an image extension
called "AREA" in the science data product. For imaging mode exposures, the
values of the PIXAR_SR and PIXAR_A2 keywords in the photom reference table
will also be copied into keywords of the same name in the primary header of
the science data product.

Note that throughout all of this process the pixel values of the science data
product are not actually changed. All of the reference data is simply attached
to the science data product in some way.

Upon successful completion of this step, the status keyword S_PHOTOM will be
set to COMPLETE.
