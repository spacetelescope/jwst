.. _readnoise_reffile:

READNOISE Reference File
------------------------

:REFTYPE: READNOISE
:Data model: `~jwst.datamodels.ReadnoiseModel`

The READNOISE reference file contains a pixel-by-pixel map of read noise,
which is used in estimating the expected noise in each pixel.

.. _readnoise_selectors:

Reference Selection Keywords for READNOISE
++++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate READNOISE references based on the following keywords.
READNOISE is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ==========================================================
Instrument Keywords
========== ==========================================================
FGS        INSTRUME, DETECTOR, READPATT, SUBARRAY, DATE-OBS, TIME-OBS
MIRI       INSTRUME, DETECTOR, READPATT, SUBARRAY, DATE-OBS, TIME-OBS
NIRCam     INSTRUME, DETECTOR, READPATT, SUBARRAY, DATE-OBS, TIME-OBS
NIRISS     INSTRUME, DETECTOR, READPATT, SUBARRAY, DATE-OBS, TIME-OBS
NIRSpec    INSTRUME, DETECTOR, READPATT, SUBARRAY, DATE-OBS, TIME-OBS
========== ==========================================================

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
The FITS primary HDU does not contain a data array.
The characteristics of the FITS extensions are as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
=======  ========  =====  ==============  =========

