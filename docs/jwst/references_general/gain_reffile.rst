.. _gain_reffile:

GAIN reference file
-------------------

:REFTYPE: GAIN
:Data model: `~jwst.datamodels.GainModel`

The GAIN reference file contains a pixel-by-pixel gain map, which can be
used to convert pixel values from units of DN to electrons. The gain
values are assumed to be in units of e/DN.

.. _gain_selectors:

Reference Selection Keywords for GAIN
+++++++++++++++++++++++++++++++++++++
CRDS selects appropriate GAIN references based on the following keywords.
GAIN is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ================================================
Instrument Keywords
========== ================================================
FGS        INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS
MIRI       INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS
NIRCam     INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS
NIRISS     INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS
NIRSpec    INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS
========== ================================================

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for GAIN
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in GAIN reference files,
because they are used as CRDS selectors
(see :ref:`gain_selectors`):

===============  ==============================
Keyword          Data Model Name
===============  ==============================
DETECTOR         model.meta.instrument.detector
SUBARRAY         model.meta.subarray.name
BUNIT\ :sup:`1`  model.meta.bunit_data
===============  ==============================

:sup:`1` BUNIT is not used as a CRDS selector, but is required in the
"SCI" extension header of GAIN reference files to document the units
of the data. The expected value is "ELECTRONS/DN".

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

