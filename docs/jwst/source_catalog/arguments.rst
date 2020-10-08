Arguments
=========

The ``source_catalog`` step uses the following user-settable arguments:

* ``--bkg_boxsize``: A floating-point value giving the background mesh box size in pixels

* ``--kernel_fwhm``: A floating-point value giving the Gaussian kernel
  FWHM in pixels [default=2.0]

* ``--snr_threshold``: A floating-point value that sets the SNR
  threshold (above background) for source detection [default=3.0]

* ``--npixels``: An integer value that sets the minimum number of
  pixels in a source [default=5]

* ``--deblend``: A boolean indicating whether to deblend sources
  [default=False]

* ``--aperture_ee1``: An integer value of the smallest aperture encircled energy value [default=30]

* ``--aperture_ee2``: An integer value of the middle aperture encircled energy value [default=50]

* ``--aperture_ee3``: An integer value of the largest aperture encircled energy value [default=70]

* ``--ci1_star_threshold``: A floating-float value of the threshold compared to the concentration
  index calculated from aperture_ee1 and aperture_ee2 that is used to determine whether a source
  is a star. Sources must meet the criteria of both ci1_star_threshold and ci2_star_threshold to be
  considered a star.

* ``--ci2_star_threshold``: A floating-float value of the threshold compared to the concentration
  index calculated from aperture_ee2 and aperture_ee3 that is used to determine whether a source
  is a star. Sources must meet the criteria of both ci1_star_threshold and ci2_star_threshold to be
  considered a star.

* ``--suffix``: A string value giving the file name suffix to use for
  the output catalog file [default='cat']
