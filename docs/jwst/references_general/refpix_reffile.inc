.. _refpix_reffile:

REFPIX Reference File
---------------------

:REFTYPE: REFPIX
:Data model: `~jwst.datamodels.IRS2Model`

The REFPIX reference file contains the complex coefficients
for the correction.

.. include:: ../references_general/refpix_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for REFPIX
+++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in REFPIX reference files,
because they are used as CRDS selectors
(see :ref:`refpix_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
READPATT   model.meta.exposure.readpatt
=========  ==============================

Reference File Format
+++++++++++++++++++++
REFPIX reference files are FITS format, with 1 BINTABLE extension.
The FITS primary HDU does not contain a data array.
The BINTABLE extension is labeled with EXTNAME = "IRS2" and has the
following column characteristics:

=======  =========
Column   Data type
=======  =========
alpha_0  float32
alpha_1  float32
alpha_2  float32
alpha_3  float32
beta_0   float32
beta_1   float32
beta_2   float32
beta_3   float32
=======  =========

The "alpha" arrays contain correction multipliers to the data,
and the "beta" arrays contain correction multipliers to the reference
output. Both arrays have 4 components; one for each sector.
The coefficients are intrinsically complex values, but have their
real and imaginary parts stored in alternating table rows, i.e. row 1
contains the real components of all coefficients and row 2 contains
the corresponding imaginary components for each.
This storage scheme results in a total of 2916352 (2048 * 712 * 2)
rows in the table.
