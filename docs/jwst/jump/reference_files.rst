Reference File Types
=====================

The `jump` step uses two reference files: GAIN and READNOISE.
Both are necessary for proper computation of noise estimates.

GAIN reference file
-------------------

:REFTYPE: GAIN
:Data model: `GainModel`

The GAIN reference file contains a pixel-by-pixel gain map, which is used
to temporarily convert pixel values in the `jump` step from units of DN to
electrons. The gain values are assumed to be in units of e/DN.

.. include:: gain_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for GAIN
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in GAIN reference files,
because they are used as CRDS selectors
(see :ref:`gain_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
SUBARRAY   model.meta.subarray.name
=========  ==============================

Reference File Format
+++++++++++++++++++++
GAIN reference files are FITS files with a single IMAGE extension.
The FITS primary data array is assumed to be empty.
The characteristics of the FITS extensions are as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
=======  ========  =====  ==============  =========

READNOISE Reference File
------------------------

:REFTYPE: READNOISE
:Data model: `ReadnoiseModel`

The READNOISE reference file contains a pixel-by-pixel map of read noise,
which is used in estimating the expected noise in each pixel.

.. include:: readnoise_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for READNOISE
++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in READNOISE reference files,
because they are used as CRDS selectors
(see :ref:`readnoise_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
READPATT   model.meta.exposure.readpatt
SUBARRAY   model.meta.subarray.name
=========  ==============================

Reference File Format
+++++++++++++++++++++
The READNOISE reference file is a FITS file with a single IMAGE extension,
which contains a 2-D floating-point array of read noise values per pixel.
The units of the read noise should be DN and should be the
CDS (Correlated Double Sampling) read noise, i.e. the effective noise between
any pair of non-destructive detector reads.
The FITS primary data array is assumed to be empty.
The characteristics of the FITS extensions are as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
=======  ========  =====  ==============  =========

