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

