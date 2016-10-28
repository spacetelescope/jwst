Reference Files
===============
The Jump step uses two reference files: GAIN and READNOISE. The gain values
are used to temporarily convert the pixel values from units of DN to
electrons. The read noise values are used as part of the noise estimate for
each pixel. Both are necessary for proper computation of noise estimates.

Gain Reference File
-------------------
The gain reference file is a FITS file with a single IMAGE extension,
with ``EXTNAME=SCI``, which contains a 2-D floating-point array of gain values
(in e/DN) per pixel. The REFTYPE value is ``GAIN``.

The gain reference files are selected from CRDS based on instrument, detector
and, where necessary, subarray.

Read Noise Reference File
-------------------------
The read noise reference file is a FITS file with a single IMAGE extension,
with ``EXTNAME=SCI``, which contains a 2-D floating-point array of read noise values
per pixel. The units of the read noise should be electrons and should be the
CDS (Correlated Double Sampling) read noise, i.e. the effective noise between
any pair of non-destructive detector reads. The REFTYPE value is
``READNOISE``.

The read noise reference files are
selected from CRDS based on instrument, detector, subarray, and readpatt.

