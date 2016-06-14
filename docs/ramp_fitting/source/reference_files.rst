Reference Files
===============

The ramp fitting step uses READNOISE and GAIN reference files.
The gain is used to temporarily convert the pixel values from units of
DN to electrons and to convert the results of ramp fitting back to DN.

CRDS Selection Criteria
-----------------------
The readnoise reference file is selected by instrument, detector and, where
necessary, subarray.

The gain reference file is selected based on instrument, detector and,
where necessary, subarray.

Reference File Format
---------------------
Read noise reference files are FITS format files having a single IMAGE
extension, labeled ``SCI``, which contains a 2-D floating-point array of read
noise values per pixel. The units of the read noise should be electrons and
should be the CDS (Correlated Double Sampling) read noise, i.e. the effective
noise between any pair of non-destructive detector reads.

The gain reference file is a FITS file with a single IMAGE extension,
labeled ``SCI``, which contains a 2-D floating-point array of gain values
(in e/DN) per pixel.
