Reference Files
===============
The ``straylight`` step uses the :ref:`MRSXARTCORR <mrsxartcorr_reffile>` reference
file, which stores vectors describing the appropriate cross-artifact convolution kernel
for each MRS band.  In Channel 1 these vectors include power in a broad Lorentzian core
plus a pair of double-Gaussian profiles.  In Channels 2 and 3 these vectors include only
power in the broad Lorentzian, while in Channel 4 there is no correction.

Optionally, the ``straylight`` step can also use the :ref:`REGIONS <regions_reffile>`
reference files during correction of residual cosmic ray showers.  These reference files
describe the footprint of the dispersed IFU traces on the detector and are used to
identify inter-slice pixels that do not see direct light from the sky.

.. include:: ../references_general/mrsxartcorr_reffile.inc
