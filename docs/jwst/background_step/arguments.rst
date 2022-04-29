.. _bkg_step_args:

Step Arguments
==============
The background image subtraction step has four optional arguments.
The first two are used only when the step is applied to non-WFSS exposures.
They are used in the process of creating an average background image, to
control the sigma clipping, and are passed as arguments to the astropy
``sigma_clip`` function:

``--sigma``
  When combining multiple background images, the number of standard deviations
  to use for the clipping limit.
  Defaults to 3.

``--maxiters``
  When combining multiple background images, the number of clipping iterations
  to perform, or ``None`` to clip until convergence is achieved.
  Defaults to ``None``.

``--save_combined_background``
  Saves the combined background image used for background subtraction.
  Defaults to ``False``.

``--wfss_mmag_extract``
  Only applies to Wide Field Slitless Spectroscopy (WFSS) exposures.
  Sets the minimum (faintest) magnitude limit to use when selecting sources
  from the WFSS source catalog, based on the value of `isophotal_abmag` in the
  source catalog. Defaults to ``None``.
