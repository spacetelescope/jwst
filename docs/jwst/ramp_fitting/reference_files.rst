Reference Files
===============
The ``ramp_fit`` step uses two reference file types: GAIN and READNOISE.
During ramp fitting, the GAIN values are used to temporarily convert the pixel
values from units of DN to electrons, and convert the results of ramp fitting
back to DN.  The READNOISE values are used as part of the noise estimate for
each pixel. Both are necessary for proper computation of noise estimates.

.. include:: ../gain_reffile/gain_reference_file.rst

.. include:: ../readnoise_reffile/readnoise_reference_file.rst

