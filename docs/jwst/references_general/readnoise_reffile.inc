READNOISE Reference File
------------------------

:REFTYPE: READNOISE
:Data model: `~jwst.datamodels.ReadnoiseModel`

The READNOISE reference file contains a pixel-by-pixel map of read noise,
which is used in estimating the expected noise in each pixel.

.. include:: ../references_general/readnoise_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for READNOISE
++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in READNOISE reference files,
because they are used as CRDS selectors
(see :ref:`readnoise_selectors`):

===============  ==============================
Keyword          Data Model Name
===============  ==============================
DETECTOR         model.meta.instrument.detector
READPATT         model.meta.exposure.readpatt
SUBARRAY         model.meta.subarray.name
BUNIT\ :sup:`1`  model.meta.bunit_data
===============  ==============================

:sup:`1` BUNIT is not used as a CRDS selector, but is required in the
"SCI" extension header of READNOISE reference files to document the units
of the data. The expected value is "DN".

Reference File Format
+++++++++++++++++++++
The READNOISE reference file is a FITS file with a single IMAGE extension,
which contains a 2-D floating-point array of read noise values per pixel.
**The units of the read noise should be DN and should be the
CDS (Correlated Double Sampling) read noise**, i.e. the effective noise between
any pair of non-destructive detector reads.
The FITS primary HDU does not contain a data array.
The characteristics of the FITS extensions are as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
=======  ========  =====  ==============  =========

