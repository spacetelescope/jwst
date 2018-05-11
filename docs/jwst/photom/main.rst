Description
============

The photom step loads - and in some cases applies - information into a
data product that allows for the
conversion of count rates to absolute flux units. The flux conversion
information is read from the photometric reference file. The exact nature
of the information that's stored in the reference file and loaded into the
science data product depends on the instrument mode.

Upon successful completion of this step, the status keyword S_PHOTOM will be
set to COMPLETE.

Imaging and non-IFU Spectroscopy
--------------------------------

Photom Data
^^^^^^^^^^^

For these instrument modes the photom reference file contains a table of
exposure parameters that define various instrument configurations and the flux
conversion data for each of those configurations. The table contains one row
for each allowed combination of exposure parameters,
such as detector, filter, pupil, and grating. The photom step searches the
table for the row that matches the parameters of the science exposure and
then copies the calibration information from that table row into the science
product. Note that for NIRSpec fixed-slit mode, the step will search the table
for each slit in use in the exposure, using the table row that corresponds to
each slit.

For these table-based reference files, the calibration information in each row
includes a scalar flux conversion constant, as well as optional arrays of
wavelength and relative response (as a function of wavelength). The scalar
conversion constant in a selected
table row is copied into the keyword PHOTMJSR in the primary header of the
science product. The value of PHOTMJSR can then be used to convert data from
units of DN/sec to MJy/steradian. The step also computes, on the fly,
the equivalent conversion factor for converting the data to units of
microJy/square-arcsecond and stores this value in the header keyword PHOTUJA2.

If the photom step finds that the wavelength and relative response arrays are
populated in the selected table row, it copies those arrays to a table extension
called "RELSENS" in the science data product.

None of the conversion factors are actually applied to the data for these
observing modes. They are simply attached to the science product.

Pixel Area Data
^^^^^^^^^^^^^^^

For imaging modes, the photom step loads data from a pixel area map
reference file and appends it to the science data product. The 2D
data array from the pixel area map is copied into an image extension
called "AREA" in the science data product.

The process of attaching the pixel
area data also populates the keywords PIXAR_SR and PIXAR_A2 in the primary
header of the science product, which give the average pixel area in units of
steradians and square arcseconds, respectively.
Both the photom and pixel area reference files contain the average pixel
area values in their primary headers. The photom step copies the values from
the pixel area reference file to populate the PIXAR_SR and PIXAR_A2 keywords
in the science data. It will issue a warning if the values of those keywords
in the two reference files differ by more than 0.1%.

NIRSpec IFU
-----------

The photom step uses the same type of tabular reference file for NIRSpec IFU
exposures as discussed above for other modes, where there is a single table
row that corresponds to a given exposure's filter and grating settings. It
retreives the scalar conversion constant, as well as the 1D wavelength and
relative response arrays, from that row. It also loads the IFU pixel area
data from the pixel area reference file.

It then uses the scalar conversion constant, the 1D wavelength and relative
response, and pixel area data to compute a 2D sensitivity map (pixel-by-pixel)
for the entire 2D science image. The 2D SCI and ERR arrays in the science
exposure are divided by the 2D sensitivity map, which converts the science
pixels from units of DN/sec to mJy/arcsec\ :sup:`2`\ . Furthermore, the
2D sensitivity array is stored in a new extension of the science exposure
called "RELSENS2D". The BUNIT keyword value in the SCI and ERR extension
headers of the science product are updated to reflect the change in units.

MIRI MRS
--------

For the MIRI MRS mode, the photom reference file contains 2D arrays of sensitivity
factors and pixel sizes that are loaded into the step. As with NIRSpec IFU, the
sensitivity and pixel size data are used to compute a 2D sensitivity map
(pixel-by-pixel) for the entire science image. This is divided into both
the SCI and ERR arrays of the science exposure, which converts the pixel values
from units of DN/sec to mJy/arcsec\ :sup:`2`\ . The 2D sensitivity array is
also stored in a "RELSENS2D" extension of the science exposure.
The BUNIT keyword value in the SCI and ERR extension
headers of the science product are updated to reflect the change in units.

