Step Arguments
==============
The ``badpix_selfcal`` step has the following optional arguments.

``--flagfrac_lower`` (float, default=0.001)
  The fraction of pixels to flag as outliers on the low-flux 
  side of the smoothed-subtracted image.

``--flagfrac_upper`` (float, default=0.001)
  The fraction of pixels to flag as outliers on the high-flux
  side of the smoothed-subtracted image.

``--kernel_size`` (integer, default=15)
  The size of the kernel to use for the median filter, which is applied 
  in the spectral direction to make the smoothed image.

``--save_flagged_bkg`` (boolean, default=False)
  Whether to save the flagged background exposures to fits files.
