.. _bkg_step_args:

Step Arguments
==============
The background image subtraction step has two optional arguments, both of
which are used only when the step is applied to non-WFSS exposures.
They are used in the process of creating an average background image, to
control the sigma clipping, and are passed as arguments to the astropy
``sigma_clip`` function:

``--sigma``
  The number of standard deviations to use for the clipping limit.
  Defaults to 3.

``--maxiters``
  The number of clipping iterations to perform, or ``None`` to clip until
  convergence is achieved. Defaults to ``None``.

``--save_combined_background``
  Saves the combined, averaged background image used for background
  subtraction. Defaults to ``False``.

``--mmag_extract``
  Only applies to Wide Field Slitless Spectroscopy (WFSS) exposures.
  Sets the minimum magnitude limit to use when selecting sources from the
  WFSS source catalog, which is used to exclude source-impacted pixels
  from the WFSS image when computing a scaling value. Defaults to 99.
