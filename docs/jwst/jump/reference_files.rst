Reference File Types
=====================

The ``jump`` step uses two reference files: GAIN and READNOISE.
The GAIN reference file is used to temporarily convert pixel values in
the ``jump`` step from units of DN to electrons.
The READNOISE reference file is used in estimating the expected noise
in each pixel.
Both are necessary for proper computation of noise estimates within the
``jump`` step.

:ref: `GAIN <../references_general/gain_reffile.rst>`

:ref: `READNOISE <../references_general/readnoise_reffile.rst>`
