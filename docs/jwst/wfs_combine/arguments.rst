.. _wfs_combine_args:

Step Arguments
==============
The ``wfs_combine`` step has the following optional arguments.

``--do_refine`` (boolean, default=False)
  If set to `True`, the nominal image offsets computed from the WCS information are
  refined using image cross-correlation. See the :ref:`wfs_comb_algo` section for details.

``--flip_dithers`` (boolean, default=True)
  When set to `True`, the output star in the combined image from the pairs of WFS images will
  always be at the same pixel location.

``--psf_size`` (float, default=100)
  The largest PSF size in pixels to use for the alignment.
  This is only used when ``do_refine`` is `True`.

``--blur_size`` (float, default=10)
  The smoothing that is applied for the initial centroiding.
  This is only used when ``do_refine`` is `True`.

``--n_size`` (int, default=2)
  This controls the size of the box used to interpolate in the input images.
  Should never need to be changed.
