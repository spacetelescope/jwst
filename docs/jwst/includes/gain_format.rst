GAIN Reference Files
--------------------

.. include:: ../includes/standard_keywords.rst

Reference Selection Keywords for GAIN
-------------------------------------
CRDS selects appropriate GAIN references based on the following keywords.
GAIN is not applicable for instruments not in the table.
Non-standard keywords used for file selection are *required*.

========== ================================================
Instrument Keywords                                         
========== ================================================
FGS        INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS 
MIRI       INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS 
NIRCAM     INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS 
NIRISS     INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS 
NIRSPEC    INSTRUME, DETECTOR, SUBARRAY, DATE-OBS, TIME-OBS 
========== ================================================

Type Specific Keywords for GAIN
-------------------------------
==========  =========================
Keyword     Data models
==========  =========================
GAINFACT    meta.exposure.gain_factor
==========  =========================

The value of GAINFACT is expected to be in the range 1..10 for NIRSPEC.

GAIN Reference Format
---------------------

The gain reference file is a FITS file with a single IMAGE extension, with
``EXTNAME=SCI``, which contains a 2-D floating-point array of gain values (in
e/DN) per pixel. The REFTYPE value is ``GAIN``.

