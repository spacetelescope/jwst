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

``--soss_source_percentile``
  The threshold flux percentile, above which values are deemed to be source or contaminated.
  Pixels with flux below this percentile will be added to the background mask. The
  default value is 35.0.

``--soss_bkg_percentile``
  This pair of percentile values describes the range of flux percentiles in the
  background mask to use for reference template scaling. The default is [25.0, 50.0].

``--wfss_mmag_extract``
  Only applies to Wide Field Slitless Spectroscopy (WFSS) exposures.
  Sets the minimum (faintest) magnitude limit to use when selecting sources
  from the WFSS source catalog, based on the value of `isophotal_abmag` in the
  source catalog. Defaults to ``None``.

``--wfss_maxiter``
  Only applies to Wide Field Slitless Spectroscopy (WFSS) exposures.
  Sets the maximum number of iterations allowed for iterative outlier rejection
  during determination of the reference background scaling factor. Defaults to 5.

``--wfss_rms_stop``
  Only applies to Wide Field Slitless Spectroscopy (WFSS) exposures.
  If the percentage difference in the RMS of the background-subtracted image
  between iterations is smaller than this value, stop the iterative outlier
  rejection process.
  Defaults to 0, i.e., do all iterations up to ``wfss_maxiter``.

``--wfss_outlier_percent``
  Only applies to Wide Field Slitless Spectroscopy (WFSS) exposures.
  Sets the percentile of outliers in the data to reject on both the low and high end
  per iteration during determination of the reference background scaling factor
  Defaults to 1, i.e., keep the middle 98 percent of the data each iteration.

``--bkg_list``
  Provides a list of background files to combine and use for subtraction. It can
  have one or more files separated by a comma with no spaces. This argument will
  be ignored for WFSS or if an association file is provided as input for the
  step.
