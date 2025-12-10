.. _picture_frame_arguments:

Step Arguments
==============

The ``picture_frame`` step has the following optional arguments to control
the behavior of the processing.

``--mask_science_regions`` (boolean, default=True)
  Mask regions of the image defined by WCS bounding
  boxes for slits/slices, as well as any regions known to be
  affected by failed-open MSA shutters.

``--n_sigma`` (float, default=2.0)
  The sigma-clipping threshold to use when searching for outliers
  and illuminated pixels to be excluded from use in the artifact
  scaling process.

``--save_mask`` (boolean, default=False)
  If set, the background mask constructed by the step will be saved to
  a file with suffix "pctfrm_mask".

``--save_correction`` (boolean, default=False)
  If set, the correction scaled to and subtracted from the input data
  will be saved to a file with suffix "pctfrm_correction".
