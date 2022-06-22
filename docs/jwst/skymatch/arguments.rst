.. _skymatch_arguments:

Step Arguments
==============
The ``skymatch`` step uses the following optional arguments:

**General sky matching parameters:**

``skymethod`` (str, default='match')
  The sky computation algorithm to be used.
  Allowed values: `local`, `global`, `match`, `global+match`

``match_down`` (boolean, default=True)
  Specifies whether the sky *differences* should be subtracted from images with
  higher sky values (``match_down`` = `True`) in order to match the image with the
  lowest sky, or sky differences should be added to the images with lower sky
  values to match the sky of the image with the highest sky value
  (``match_down`` = `False`). **NOTE**: this argument only applies when
  ``skymethod`` is either `match` or `global+match`.

``subtract`` (boolean, default=False)
  Specifies whether the computed sky background values are to be subtracted from
  the images. The BKGSUB keyword (boolean) will be set in each output image to
  record whether or not the background was subtracted.

**Image bounding polygon parameters:**

``stepsize`` (int, default=None)
  Spacing between vertices of the images bounding polygon. The default value of
  `None` creates bounding polygons with four vertices corresponding to the corners
  of the image.

**Sky statistics parameters:**

``skystat`` (str, default='mode')
  Statistic to be used for sky background
  computations. Supported values are: `mean`, `mode`, `midpt`,
  and `median`.

``dqbits`` (str, default='~DO_NOT_USE+NON_SCIENCE')
  The DQ bit values from the input images' DQ arrays that
  should be considered "good" when building masks for sky computations. See
  DQ flag :ref:`dq_parameter_specification` for details. The default value
  rejects pixels flagged as either 'DO_NOT_USE' or 'NON_SCIENCE' and considers
  all others to be good.

``lower`` (float, default=None)
  An optional value indicating the lower limit of usable pixel
  values for computing the sky. This value should be specified in the units
  of the input images.

``upper`` (float, default=None)
  An optional value indicating the upper limit of usable pixel
  values for computing the sky. This value should be specified in the units
  of the input images.

``nclip`` (int, default=5)
  The number of clipping iterations to use when computing sky values.

``lsig`` (float, default=4.0)
  Lower clipping limit, in sigma, used when computing the sky value.

``usig`` (float, default=4.0)
  Upper clipping limit, in sigma, used when computing the sky value.

``binwidth`` (float, default=0.1)
  Bin width, in sigma, used to sample the distribution of pixel
  values in order to compute the sky background using statistics
  that require binning, such as `mode` and `midpt`.
