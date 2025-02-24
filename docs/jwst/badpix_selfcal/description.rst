Description
===========

:Class: `jwst.badpix_selfcal.BadpixSelfcalStep`
:Alias: badpix_selfcal

Overview
--------
The ``badpix_selfcal`` step flags bad pixels in the input data using a self-calibration 
technique based on median filtering along the spectral axis. 
When additional exposures are available, those are used in combination with the science
exposure to identify bad pixels; when unavailable, the step will be skipped with a warning
unless the ``force_single`` parameter is set True. In that case, the science data alone is
used as its own "background".
This correction is applied to `~jwst.datamodels.IFUImageModel` data
directly after the :ref:`assign_wcs <assign_wcs_step>` correction has been applied
in the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline.

Input details
-------------
The input data must be in the form of a `~jwst.datamodels.IFUImageModel` or
a `~jwst.datamodels.ModelContainer` containing exactly one
science exposure and any number of additional exposures.
A fits or association file 
that can be read into one of these data models is also acceptable.
Any exposure with the metadata attribute ``asn.exptype`` set to 
``background`` or ``selfcal`` will be used in conjunction with the science
exposure to construct the combined background image. 

Algorithm
---------
The algorithm relies on the assumption that bad pixels are outliers in the data along
the spectral axis. The algorithm proceeds as follows:

* A combined background image is created. If additional (``selfcal`` or ``background``)
  exposures are available, 
  the pixelwise minimum of all background, selfcal, and science exposures is taken. 
  If no additional exposures are available, the science data itself is passed in 
  without modification, serving as the "background image" for the rest of the procedure, 
  i.e., true self-calibration.
* For MIRI MRS, any residual pedestal in the effective dark current is subtracted from
  each exposure using the unilluminated region of the detector between spectral channels
  prior to performing the pixelwise minimum combination.
* The combined background image is median-filtered, ignoring NaNs, along the spectral axis 
  with a user-specified kernel size. The default kernel size is 15 pixels.
* The difference between the original background image and the median-filtered background image
  is taken. The highest- and lowest-flux pixels in this difference image are
  flagged as bad pixels. The default fraction of pixels to flag is 0.1% of the total number of pixels
  on each of the high-flux and low-flux ends of the distribution. These fractions can be adjusted
  using the ``flagfrac_lower`` and ``flagfrac_upper`` parameters for the low- and high-flux ends
  of the distribution, respectively. The total fraction of flagged pixels is thus 
  ``flagfrac_lower + flagfrac_upper``.
* The bad pixels are flagged in the input data by setting the DQ flag to
  "OTHER_BAD_PIXEL" and "DO_NOT_USE".
* The bad pixels are also flagged in each exposure with ``asn.exptype`` equal to ``background``,
  if available.

Output product
--------------
The output is a new copy of the input `~jwst.datamodels.IFUImageModel`, with the
bad pixels flagged.  If the entire ``calwebb_spec2`` pipeline is run, the background
exposures passed into the ``background`` step will include the flags from the
``badpix_selfcal`` step.
