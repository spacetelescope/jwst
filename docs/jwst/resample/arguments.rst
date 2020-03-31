.. _resample_step_args:

Step Arguments
==============
The `resample` step has the following optional arguments that control
the behavior of the processing and the characteristics of the resampled
image.

``--pixfrac`` (float, default=1.0)
  The fraction by which input pixels are "shrunk" before being drizzled
  onto the output image grid, given as a real number between 0 and 1.

``--kernel`` (str, default='square')
  The form of the kernel function used to distribute flux onto the output
  image.

``--fillval`` (str, default='INDEF')
  The value to assign to output pixels that have zero weight or do not
  receive any flux from any input pixels during drizzling.

``--weight_type`` (str, default='exptime')
  The weighting factor for each input image. If `weight_type=exptime`,
  the scaling value will be set equal to the exposure time found in
  the image header.

``--single`` (bool, default=False)
  Resample each input image into a separate output.

``--blendheaders`` (bool, default=True)
  Apply `blendmodels` on all of the input images to combine ('blend')
  their meta data into the output resampled image.
