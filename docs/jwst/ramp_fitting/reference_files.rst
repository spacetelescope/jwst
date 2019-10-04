Reference Files
===============
The ``ramp_fit`` step uses two reference file types: :ref:`GAIN <gain_reffile>`
and :ref:`READNOISE <readnoise_reffile>`.
During ramp fitting, the GAIN values are used to temporarily convert the pixel
values from units of DN to electrons, and convert the results of ramp fitting
back to DN.  The READNOISE values are used as part of the noise estimate for
each pixel. Both are necessary for proper computation of noise estimates.

:ref:`GAIN <gain_reffile>`

:ref:`READNOISE <readnoise_reffile>`

