Description
-----------
The ``hlsp`` step is one of the coronagraphic-specific steps in the ``coron``
sub-package. It produces high-level science products for KLIP-processed
(PSF-subtracted) coronagraphic images. The step is currently a prototype and
produces two simple products: a signal-to-noise ratio (SNR) image and a table
of contrast data. The SNR image is computed by simply taking the ratio of the
SCI and ERR arrays of the input target image. The contrast
data are in the form of azimuthally-averaged noise versus radius. The noise
is computed as the 1-sigma standard deviation within a set of concentric
annuli centered in the input image. The annuli regions are computed to the
nearest whole pixel; no sub-pixel calculations are performed.

.. Note:: This step is not currently included in the :ref:`calwebb_coron3 <calwebb_coron3>`
   pipeline, but can be run standalone.

Arguments
---------
The ``hlsp`` step has one optional argument:

``--annuli_width``
  which is an integer parameter with a default value of 2 and is used to
  specify the width, in pixels, of the annuli to use when computing the contrast
  curve data.

Inputs
------

2D image
^^^^^^^^
:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _psfsub

The input is the KLIP-processed (PSF-subtracted) image to be analyzed.

Outputs
-------

2D SNR image
^^^^^^^^^^^^
:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _snr

The computed SNR image.

Contrast table
^^^^^^^^^^^^^^
:Data model: `~jwst.datamodels.ContrastModel`
:File suffix: _contrast

The table of contrast data, containing columns of radii (in pixels) and 1-sigma noise.

Reference Files
---------------
The ``hlsp`` step does not use any reference files.
