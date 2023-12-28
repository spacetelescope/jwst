.. _nsclean_arguments:

Step Arguments
==============

The ``nsclean`` step has the following optional arguments to control
the behavior of the processing.

``--mask_spectral_regions`` (boolean, default=True)
  Mask regions in IFU and MOS images that are within the bounding boxes
  for each slice or slitlet defined in the WCS object of the image.

``--n_sigma`` (float, default=5.0)
  The sigma-clipping threshold to use when searching for outliers
  and illuminated pixels to be excluded from use in the fitting
  process.

``--save_mask`` (boolean, default=False)
  A flag to indicate whether the mask constructed by the step
  should be saved to a file.

``--user_mask`` (string, default=None)
  Path to a user-supplied mask file. If supplied, the mask is used
  directly and the process of creating a mask in the step is skipped.
