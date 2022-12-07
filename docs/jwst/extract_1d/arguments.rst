Step Arguments
==============

The ``extract_1d`` step has the following step-specific arguments.

``--smoothing_length``
  If ``smoothing_length`` is greater than 1 (and is an odd integer), the
  image data used to perform background extraction will be smoothed in the
  dispersion direction with a boxcar of this width.  If ``smoothing_length``
  is None (the default), the step will attempt to read the value from the
  EXTRACT1D reference file.  If a value is specified in the reference file,
  that value will be used.  Note that by specifying this parameter in the
  EXTRACT1D reference file a different value can be designated for each slit.
  If no value is specified either by the user or in the EXTRACT1D reference
  file, no background smoothing is done.

``--bkg_fit``
  The type of fit to perform to the background data in each image column
  (or row, if the dispersion is vertical). There are four allowed values:
  "poly", "mean", and "median", and None (the default value). If left as None,
  the step will search the reference file for a value - if none is found,
  ``bkg_fit`` will be set to "poly". If set to "poly", the background
  values for each pixel within all background regions in a given column (or
  row) will be fit with a polynomial of order "bkg_order" (see below).
  Values of "mean" and "median" compute the simple average and median,
  respectively, of the background region values in each column (or row).
  This parameter can also be specified in the EXTRACT1D reference file. If
  specified in the reference file and given as an argument to the step by
  the user, the user-supplied value takes precedence.

``--bkg_order``
  The order of a polynomial function to be fit to the background
  regions.  The fit is done independently for each column (or row, if the
  dispersion is vertical) of the input image, and the fitted curve will be
  subtracted from the target data.  ``bkg_order`` = 0 (the minimum allowed
  value) means to fit a constant.  The user-supplied value (if any)
  overrides the value in the EXTRACT1D reference file.  If neither is specified, a
  value of 0 will be used. If a sufficient number of valid data points is
  unavailable to construct the polynomial fit, the fit will be forced to
  0 for that particular column (or row). If "bkg_fit" is not "poly", this
  parameter will be ignored.

``--bkg_sigma_clip``
  The background values will be sigma-clipped to remove outlier values from
  the determination of the background. The default value is a 3.0 sigma clip.

``--log_increment``
  Most log messages are suppressed while looping over integrations, i.e. when
  the input is a CubeModel or a 3-D SlitModel.  Messages will be logged while
  processing the first integration, but since they would be the same for
  every integration, most messages will only be written once.  However, since
  there can be hundreds or thousands of integrations, which can take a long
  time to process, it would be useful to log a message every now and then to
  let the user know that the step is still running.

  ``log_increment`` is an integer, with default value 50.  If it is greater
  than 0, an INFO message will be printed every ``log_increment``
  integrations, e.g. "... 150 integrations done".

``--subtract_background``
  This is a boolean flag to specify whether the background should be
  subtracted.  If None, the value in the EXTRACT1D reference file (if any)
  will be used.  If not None, this parameter overrides the value in the
  reference file.

``--use_source_posn``
  This is a boolean flag to specify whether the target and background extraction
  region locations specified in the EXTRACT1D reference file should be shifted
  to account for the expected position of the source. If None (the default),
  the step will make the decision of whether to use the source position based
  on the observing mode and the source type. The source position will only be
  used for point sources and for modes where the source could be located
  off-center due to things like nodding or dithering. If turned on, the sky
  (RA/Dec) position of the source is used in conjunction with the World
  Coordinate System (WCS) to compute the x/y source location. For long-slit
  type modes (e.g. MIRI LRS and NIRSpec fixed-slit and MOS), only the position
  in the cross-dispersion direction is used to potentially offset the
  extraction regions in that direction.

``--center_xy``
  A list of two integer values giving the desired x/y location for the center
  of the circular extraction aperture used for extracting spectra from 3-D
  IFU cubes. Ignored for non-IFU modes and non-point sources. Must be given in
  x,y order and in units of pixels along the x,y axes of the 3-D IFU cube, e.g.
  ``--center_xy="27,28"``. If given, the values override any position derived
  from the use of the ``use_source_posn`` argument. Default is None.

``--apply_apcorr``
  Switch to select whether or not to apply an APERTURE correction during the
  Extract1dStep processing. Default is ``True``

``--soss_atoca``
  This is a NIRISS-SOSS algorithm-specific parameter; if True, use the ATOCA
  algorithm to treat order contamination. Default is ``True``.

``--soss_threshold``
  This is a NIRISS-SOSS algorithm-specific parameter; this sets the threshold
  value for a pixel to be included when modelling the spectral trace. The default
  value is 0.01.

``--soss_n_os``
  This is a NIRISS-SOSS algorithm-specific parameter; this is an integer that sets
  the oversampling factor of the underlying wavelength grid used when modeling the
  trace. The default value is 2.

``--soss_estimate``
  This is a NIRISS-SOSS algorithm-specific parameter; filename or SpecModel of the
  estimate of the target flux. The estimate must be a SpecModel with wavelength and
  flux values.

``--soss_wave_grid_in``
  This is a NIRISS-SOSS algorithm-specific parameter; filename or SossWaveGridModel
  containing the wavelength grid used by ATOCA to model each valid pixel of the
  detector. If not given, the grid is determined based on an estimate of the flux
  (soss_estimate), the relative tolerance (soss_rtol) required on each pixel model
  and the maximum grid size (soss_max_grid_size).

``--soss_wave_grid_out``
  This is a NIRISS-SOSS algorithm-specific parameter; filename to hold the wavelength
  grid calculated by ATOCA, stored in a SossWaveGridModel.

``--soss_rtol``
  This is a NIRISS-SOSS algorithm-specific parameter; the relative tolerance needed on a
  pixel model. It is used to determine the sampling of the soss_wave_grid when not
  directly given. Default value is 1.e-4.

``--soss_max_grid_size``
  This is a NIRISS-SOSS algorithm-specific parameter; the maximum grid size allowed. It is
  used when soss_wave_grid is not provided to make sure the computation time or the memory
  used stays reasonable. Default value is 20000.

``--soss_transform``
  This is a NIRISS-SOSS algorithm-specific parameter; this defines a rotation to
  apply to the reference files to match the observation. It should be specified as
  a list of three floats, with default values of None.

``--soss_tikfac``
  This is a NIRISS-SOSS algorithm-specific parameter; this is the regularization
  factor used in the SOSS extraction. If not specified, ATOCA will calculate a
  best-fit value for the Tikhonov factor.

``--soss_width``
  This is a NIRISS-SOSS algorithm-specific parameter; this specifies the aperture
  width used to extract the 1D spectrum from the decontaminated trace. The default
  value is 40.0 pixels.

``--soss_bad_pix``
  This is a NIRISS-SOSS algorithm-specific parameter; this parameter sets the method
  used to handle bad pixels. There are currently two options: "model" will replace
  the bad pixel values with a modeled value, while "masking" will omit those pixels
  from the spectrum. The default value is "model".

``--soss_modelname``
  This is a NIRISS-SOSS algorithm-specific parameter; if set, this will provide
  the optional ATOCA model output of traces and pixel weights, with the filename
  set by this parameter. By default this is set to None and this output is
  not provided.
