.. _psfmask_reffile:

PSFMASK Reference File
^^^^^^^^^^^^^^^^^^^^^^

:REFTYPE: PSFMASK
:Data model: `~jwst.datamodels.PsfMaskModel`

The PSFMASK reference file contains a 2-D mask that's used as a weight
function when computing shifts between images.

.. include:: ../references_general/psfmask_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for PSFMASK
+++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in PSFMASK reference files,
because they are used as CRDS selectors
(see :ref:`psfmask_selectors`):

=========  =================================
Keyword    Data Model Name
=========  =================================
FILTER     model.meta.instrument.filter
CORONMSK   model.meta.instrument.coronagraph
SUBARRAY   model.meta.subarray.name
=========  =================================

Reference File Format
+++++++++++++++++++++
PSFMASK reference files are FITS format, with 1 IMAGE extension.
The FITS primary HDU does not contain a data array.
The format and content of the file is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
=======  ========  =====  ==============  =========

The values in the ``SCI`` array give the mask values to be applied to
the images when computing relative shifts. The mask acts as a weighting
function when performing Fourier fits. The values range from zero
(full weighting) to one (pixel completely masked out).
