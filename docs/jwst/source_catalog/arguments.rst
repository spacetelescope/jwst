Arguments
=========

The ``source_catalog`` step uses the following user-settable arguments:

* ``--kernel_fwhm``: A floating-point value giving the Gaussian kernel
  FWHM in pixel [default=2.0]

* ``--kernel_xsize``: A floating-point value giving the kernel x size
  in pixels [default=5]

* ``--kernel_ysize``: A floating-point value giving the kernel y size
  in pixels [default=5]

* ``--snr_threshold``: A floating-point value that sets the SNR
  threshold (above background) for source detection [default=3.0]

* ``--npixels``: A floating-point value that sets the minimum number of
  pixels in a source [default=5]

* ``--deblend``: A boolean indicating whether to deblend sources
  [default=False]
