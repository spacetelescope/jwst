.. _sirskernel_reffile:

SIRSKERNEL Reference File
-------------------------

:REFTYPE: SIRSKERNEL
:Data model: `~jwst.datamodels.SIRSKernelModel`

The SIRSKERNEL reference file contains Fourier coefficients for an
optimized convolution kernel, used for correcting 1/f noise from
detector reference pixels.

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
SIRSKERNEL reference files are in ASDF format.  They contain a dictionary with
keys corresponding to lower-case detector names.  Each detector entry contains
two arrays, `gamma` and `zeta`, which contain the complex Fourier coefficients
required to apply the Simple Improved Reference Subtraction (SIRS) algorithm
for reference pixel correction.

The first dimension of the array holds the sector index (amplifier output,
4 per detector).  The second dimension matches the height of the sector, plus one.
The data type for both arrays is `complex64`.
