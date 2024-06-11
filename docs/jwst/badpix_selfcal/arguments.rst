Step Arguments
==============
The ``badpix_selfcal`` step has the following optional arguments.

``--flagfrac`` (float, default=0.001)
  The fraction of pixels to flag as outliers on each of the low-flux and high-flux
  sides of the smoothed-subtracted image.

``--kernel_size`` (integer, default=15)
  The size of the kernel to use for the median filter, which is applied 
  in the spectral direction to make the smoothed image.

``--save_flagged_bkgd`` (boolean, default=False)
  Whether to save the flagged background exposures to fits files.
