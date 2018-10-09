Reference Files
===============
The Ramp Fitting step uses two reference file types: GAIN and READNOISE.  The
GAIN reference file is also used by gain_step.

During ramp fitting, the gain values are used to temporarily convert the pixel
values from units of DN to electrons, and convert the results of ramp fitting
back to DN.  The read noise values are used as part of the noise estimate for
each pixel. Both are necessary for proper computation of noise estimates.

.. include:: ../includes/gain_format.rst

READNOISE Reference Files
-------------------------

.. include::  ../includes/standard_keywords.rst

.. include:: readnoise_selection.rst

READNOISE Format
----------------

The read noise reference file is a FITS file with a single IMAGE extension,
with ``EXTNAME=SCI``, which contains a 2-D floating-point array of read noise
values per pixel. The units of the read noise should be DN and should be the
CDS (Correlated Double Sampling) read noise, i.e. the effective noise between
any pair of non-destructive detector reads. The REFTYPE value is ``READNOISE``.


