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
REFPIX reference files are FITS format, with two BINTABLE extensions.
The FITS primary HDU does not contain a data array.
The first BINTABLE extension is labeled with EXTNAME = "IRS2" and has the
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

The "alpha" arrays contain correction multipliers to the reference output,
and the "beta" arrays contain correction multipliers to the interleaved
reference pixels. Both arrays have 4 components, one for each sector
(amplifier output).  These are applied to (i.e. multiplied by) the Fourier
transform of the interleaved reference pixel data.
The coefficients are intrinsically complex values, but have their
real and imaginary parts stored in alternating table rows, i.e. row 1
contains the real components of all coefficients and row 2 contains
the corresponding imaginary components for each.
This storage scheme results in a total of 2916352 (2048 * 712 * 2)
rows in the table.

The second BINTABLE extension is labeled with EXTNAME = "DQ" and has the
following column characteristics:

========  =========
Column    Data type
========  =========
output    int16
odd_even  int16
mask      uint32
========  =========

This table has eight rows.  The "output" column contains the amplifier output
numbers:  1, 1, 2, 2, 3, 3, 4, 4.  The "odd_even" column contains values
1 or 2, indicating that either the first or second pair of reference pixel
reads respectively should be regarded as bad.  The "mask" column contains
32-bit unsigned integer values.  The interpretation of these values was
described in the ESA CDP3 document as follows:

"There is also a DQ extension that holds a binary table with three
columns (OUTPUT, ODD_EVEN, and MASK) and eight rows. In the current
IRS2 implementation, one jumps 32 times to odd and 32 times to even
reference pixels, which are then read twice consecutively. Therefore,
the masks are 32 bit unsigned integers that encode bad interleaved
reference pixels/columns from left to right (increasing column index)
in the native detector frame. When a bit is set, the corresponding
reference data should not be used for the correction."

We assume that "native detector frame" in the above description referred to
the order that the data and interleaved reference pixels were read out from
the detector, not the physical locations of the pixels on the detector.  The
difference is that the readout direction changes when going from one
amplifier output to the next; that is, the pixels are read out from left to
right for the first and third outputs, and they are read out from right to
left for the second and fourth outputs.  Furthermore, we assume that for the
first amplifier output, it is the least significant bit in the value from the
MASK column that corresponds to the first set of four reads of interleaved
reference pixel values (reading pixels from left to right).
