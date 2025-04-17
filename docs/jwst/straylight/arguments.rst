Step Arguments
==============
The ``straylight`` step has the following optional arguments.

``--clean_showers`` (boolean, default=False)
  A flag to indicate whether to remove straylight due to residual cosmic
  ray showers using pixels between the IFU slices.

``--save_shower_model`` (boolean, default=False)
  If set, the model of the residual cosmic ray shower artifacts will
  be saved with the suffix `shower_model`.

``--shower_plane`` (int, default=3)
  Identifies the throughput plane to use from the MRS regions reference
  files to identify between-slice pixels.  Lower values identify pixels
  with less science signal, but result in fewer (or no) pixels between
  slices in some bands.

``--shower_x_stddev`` (float, default=18)
   Standard deviation to use for gaussian kernel convolution in the X
   (between-slice) direction.  Must be large enough to reach the middle
   of the IFU slices.

``--shower_y_stddev`` (float, default=5)
   Standard deviation to use for gaussian kernel convolution in the Y
   (along-slice) direction.

``--shower_low_reject`` (float, default=0.1)
   Percentile of low-valued pixels to reject when performing gaussian
   kernel convolution.  Important to ensure unmasked bad pixels do
   not bias the kernel convolution.

``--shower_high_reject`` (float, default=99.9)
   Percentile of high-valued pixels to reject when performing gaussian
   kernel convolution.  Important to ensure unmasked bad pixels do
   not bias the kernel convolution.
