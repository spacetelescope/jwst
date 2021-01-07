Reference Files
===============
The ``flat_field`` step uses four different types of reference files, depending on the
type of data being processed. Most cases just use the FLAT reference file, while NIRSpec
spectroscopic exposures use the three reference files FFLAT (fore optics),
SFLAT (spectrograph optics), and DFLAT (detector).

.. include:: ../references_general/flat_reffile.inc

Reference Files for NIRSpec Spectroscopy
===================================================
For NIRSpec spectroscopic data, the flat-field reference files allow for variations in
the flat field with wavelength, as well as from pixel to pixel.  There is a
separate flat-field reference file for each of three sections of the
instrument:  the fore optics (FFLAT), the spectrograph (SFLAT), and the
detector (DFLAT).  The contents of the reference files differ from one mode
to another (see below), but in general they may contain a flat-field image and
a 1-D array.  The image provides pixel-to-pixel values for the flat field
that may vary slowly (or not at all) with wavelength, while the 1-D array
is for a pixel-independent fast variation with wavelength. Details of the
file formats are given in the following sections.

If there is no significant slow variation with wavelength, the image will be a 2-D array;
otherwise, the image will be a 3-D array, with each plane corresponding to
a different wavelength. In the latter case, the wavelength for each plane
will be given in a table extension called WAVELENGTH in the flat-field
reference file.  The fast variation is given in a table extension called
FAST_VARIATION, with column names "slit_name", "nelem", "wavelength", and
"data" (an array of wavelength-dependent flat-field values).  Each row of
the table contains a slit name (for fixed-slit data, otherwise "ANY"), an
array of flat-field values, an array of the corresponding wavelengths, and
the number of elements ("nelem") of "data" and "wavelength" that are
populated, because the allocated array size can be larger than needed.
For some reference files there will not be any image array, in which case
all the flat field information will be taken from the FAST_VARIATION table.

The SCI extension of the reference files may contain NaNs.  If so, the
``flat_field`` step will replace these values with 1 and will flag the
corresponding pixel in the DQ extension with NO_FLAT_FIELD. The WAVELENGTH
extension is not expected to contain NaNs.

For the detector section, there is only one flat-field reference file for
each detector.  For the fore optics and the spectrograph sections, however,
there are different flat fields for fixed-slit data, IFU data, and for
multi-object spectroscopic data.  Here is a summary of the contents of these
files.

For the fore optics (FFLAT), the flat field for fixed-slit data contains just a
FAST_VARIATION table (i.e. there is no image).  This table has five rows,
one for each of the fixed slits.  The FFLAT for IFU data also contains
just a FAST_VARIATION table, but it has only one row with the value "ANY"
in the "slit_name" column.  For multi-object spectroscopic data, the FFLAT
contains four sets of images (one for each MSA quadrant), WAVELENGTH
tables, and FAST_VARIATION tables.  The images are unique to the FFLATs,
however.  The image "pixels" correspond to micro-shutter array slits, rather
than to detector pixels.  The array size is 365 columns by 171 rows,
and there are multiple planes to handle the slow variation
of flat field with wavelength.

For the spectrograph optics (SFLAT), the flat-field files have nearly the same
format for fixed-slit data, IFU, and multi-object data.  The difference is
that for fixed-slit and IFU data, the image is just a single plane,
i.e. the only variation with wavelength is in the FAST_VARIATION table,
while there are multiple planes in the image for multi-object spectroscopic
data (and therefore there is also a corresponding WAVELENGTH table, with
one row for each plane of the image).

For the detector section, the DFLAT file contains a 3-D image
(i.e. the flat field at multiple wavelengths), a corresponding
WAVELENGTH table, and a FAST_VARIATION table with one row.

As just described, there are 3 types of reference files for NIRSpec (FFLAT,
SFLAT, and DFLAT), and within each of these types, there are several formats,
which are now described.

.. include:: ../references_general/fflat_reffile.inc

.. include:: ../references_general/sflat_reffile.inc

.. include:: ../references_general/dflat_reffile.inc

