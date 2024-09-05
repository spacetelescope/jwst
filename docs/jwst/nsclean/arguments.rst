.. _nsclean_arguments:

Step Arguments
==============

The ``nsclean`` step has the following optional arguments to control
the behavior of the processing.

``--fit_method`` (str, default='fft')
  The noise fitting algorithm to use.  Options are 'fft' and 'median'.

``--fit_by_channel`` (boolean, default=False)
  If set, flicker noise is fit independently for each detector channel.
  Ignored for subarray data and for `fit_method` = 'fft'.

``--background_method`` (str, default=None)
  If 'median', the preliminary background to remove and restore
  after fitting the noise is a simple median of the background data.
  If 'model', the background data is fit with a low-resolution model
  via `~photutils.background.Background2D`.
  If None, the background value is set to 0.0.

``--background_box_size`` (list of int, default=(32,32))
  Box size for the data grid used by `~photutils.background.Background2D`
  when `background_method` = 'model'. For best results, use a
  box size that evenly divides the input image shape.

``--mask_spectral_regions`` (boolean, default=True)
  Mask regions of the image defined by WCS bounding
  boxes for slits/slices, as well as any regions known to be
  affected by failed-open MSA shutters.

``--n_sigma`` (float, default=5.0)
  The sigma-clipping threshold to use when searching for outliers
  and illuminated pixels to be excluded from use in the background
  and noise fitting processes.

``--fit_histogram`` (boolean, default=False)
  If set, the 'sigma' used with `n_sigma` for clipping outliers
  is derived from a Gaussian fit to a histogram of values.
  Otherwise, a simple iterative sigma clipping is performed.

``--single_mask`` (boolean, default=False)
  If set, a single mask will be created, regardless of
  the number of input integrations. Otherwise, the mask will
  be a 3D cube, with one plane for each integration.

``--user_mask`` (string, default=None)
  Path to a user-supplied mask file. If supplied, the mask is used
  directly and the process of creating a scene mask in the step is
  skipped.

  The mask file must contain either a `~jwst.datamodels.ImageModel`
  or a `~jwst.datamodels.CubeModel`, with image dimensions matching
  the input science data.  If an ImageModel is provided, the same
  mask will be used for all integrations.  If a CubeModel is provided,
  the number of slices must equal the number of integrations in
  the input science data.

``--save_mask`` (boolean, default=False)
  If set, the mask constructed by the step will be saved to a file
  with suffix 'mask'.

``--save_background`` (boolean, default=False)
  If set, the background fit to the group diff images will be saved
  to a file with suffix 'flicker_bkg'.

``--save_noise`` (boolean, default=False)
  If set, the residual noise fit and removed from the input data
  will be saved to a file with suffix 'flicker_noise'.
