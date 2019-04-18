Description
============

The photom step applies flux (photometric) calibrations to a data product
to convert the data from units of countrate to surface brightness.
The calibration information is read from a photometric reference file.
The exact nature of the calibration information loaded from the reference file
and applied to the science data depends on the instrument mode.

This step relies on having wavelength information available when working on
spectroscopic data (see below) and therefore the
:ref:`assign_wcs <assign_wcs_step>` step *must* be applied before executing
the photom step.

Upon successful completion of this step, the status keyword S_PHOTOM will be
set to COMPLETE.
Furthermore, the BUNIT keyword value in the SCI and ERR extension
headers of the science product are updated to reflect the change in units.

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
then loads the calibration information from that table.
Note that for NIRSpec fixed-slit mode, the step will search the table
for each slit in use in the exposure, using the table row that corresponds to
each slit.

For these table-based photom reference files, the calibration information in each
row includes a scalar flux conversion constant, as well as optional arrays of
wavelength and relative response (as a function of wavelength).
For spectroscopic data, if the photom step finds that the wavelength and relative
response arrays in the reference table row are populated, it loads those 1-D arrays
and interpolates the response values into the 2-D space of the science image based
on the wavelength at each pixel.

The combination of the scalar conversion factor and the 2-D response values are
then applied to the science data, including the SCI and ERR arrays, as well as
the variance (VAR_POISSON and VAR_RNOISE) arrays.

The scalar conversion constant is copied to the header keyword PHOTMJSR, which
gives the conversion from DN/s to MegaJy/steradian that was applied to the data.
The step also computes the equivalent conversion factor to units of
microJy/square-arcsecond and stores it in the header keyword PHOTUJA2.

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
retreives the scalar conversion constant, as well as the 1-D wavelength and
relative response arrays, from that row. It also loads the IFU pixel area
data from the pixel area reference file.

It then uses the scalar conversion constant, the 1-D wavelength and relative
response, and pixel area data to compute a 2-D sensitivity map (pixel-by-pixel)
for the entire science image. The 2-D SCI and ERR arrays in the science
exposure are divided by the 2D sensitivity map, which converts the science
pixels from units of DN/sec to mJy/arcsec\ :sup:`2`\ .
Variance arrays are divided by the square of the conversion factors.

MIRI MRS
--------

For the MIRI MRS mode, the photom reference file contains 2-D arrays of sensitivity
factors and pixel sizes that are loaded into the step. As with NIRSpec IFU, the
sensitivity and pixel size data are used to compute a 2-D sensitivity map
(pixel-by-pixel) for the entire science image. This is divided into both
the SCI and ERR arrays of the science exposure, which converts the pixel values
from units of DN/sec to mJy/arcsec\ :sup:`2`\ .
Variance arrays are divided by the square of the conversion factors.
