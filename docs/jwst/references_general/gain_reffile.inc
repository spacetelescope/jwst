GAIN reference file
-------------------

:REFTYPE: GAIN
:Data model: `~jwst.datamodels.GainModel`

The GAIN reference file contains a pixel-by-pixel gain map, which can be
used to convert pixel values from units of DN to electrons. The gain
values are assumed to be in units of e/DN.

.. include:: ../references_general/gain_selection.inc

.. include:: ../includes/standard_keywords.inc

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

