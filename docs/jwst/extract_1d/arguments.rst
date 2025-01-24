Step Arguments
==============

The ``extract_1d`` step has the following step-specific arguments.

General Step Arguments
----------------------
The following arguments apply to all modes unless otherwise specified.

``--subtract_background``
  Flag to specify whether the background should be subtracted.  If None or True,
  background subtraction will be performed if there are background regions
  specified in the reference file.  If False, no background subtraction will be
  performed.  Has no effect for NIRISS SOSS data.

``--apply_apcorr``
  Switch to select whether or not to apply an APERTURE correction during the
  Extract1dStep processing. Default is ``True``. Has no effect for NIRISS SOSS data
  or for optimal extractions.

Step Arguments for Slit and Slitless Spectroscopic Data
-------------------------------------------------------

``--extraction_type``
  Specify the extraction type.
  If 'box', a standard extraction is performed, summing over an aperture box.
  If 'optimal', a PSF-based optimal extraction is performed.
  If None, optimal extraction is attempted whenever use_source_posn is True.
  Box extraction is suitable for any input data (point sources and extended sources;
  resampled and unresampled images).  Optimal extraction is best suited for unresampled
  point sources. Currently, optimal extraction is only available for MIRI LRS Fixed Slit data.
  The default extraction type is 'box'.

``--use_source_posn``
  Specify whether the target and background extraction
  region locations specified in the :ref:`EXTRACT1D <extract1d_reffile>` reference
  file should be shifted to account for the expected position of the source. If None (the default),
  the step will decide whether to use the source position based
  on the observing mode and the source type. By default, source position corrections
  are attempted only for point sources in NIRSpec MOS/FS/BOTS and MIRI LRS fixed-slit exposures.
  Set to False to ignore position estimates for all modes; set to True to additionally attempt
  source position correction for extended sources.

``--position_offset``
  Specify a number of pixels (fractional pixels are allowed) to offset the 
  extraction aperture from the nominal position.  The default is 0.

``--model_nod_pair``
  Flag to enable fitting a negative trace during optimal extraction.
  If True, and the extraction type is 'optimal', then a negative trace
  from nod subtraction is modeled alongside the positive source during
  extraction.  This will be attempted only if the input data has been background
  subtracted and the dither pattern type indicates that 2 nods were used.
  The default value is True.

``--optimize_psf_location``
  Flag to enable PSF location optimization during optimal extraction.
  If True, and the extraction type is 'optimal', then the placement of
  the PSF model for the source location (and negative nod, if present)
  will be iteratively optimized. This parameter is recommended if
  negative nods are modeled.  The default value is True.

``--smoothing_length``
  If ``smoothing_length`` is greater than 1 (and is an odd integer), the
  image data used to perform background extraction will be smoothed in the
  dispersion direction with a boxcar of this width.  If ``smoothing_length``
  is None (the default), the step will attempt to read the value from the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file.
  If a value is specified in the reference file,
  that value will be used.  Note that by specifying this parameter in the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file, a different value can
  be designated for each slit, if desired.
  If no value is specified either by the user or in the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file,
  no background smoothing is done.

``--bkg_fit``
  The type of fit to perform to the background data in each image column
  (or row, if the dispersion is vertical). There are four allowed values:
  "poly", "mean", "median", and None (the default value). If left as None,
  the step will search the reference file for a value - if none is found,
  ``bkg_fit`` will be set to "poly". If set to "poly", the background
  values for each pixel within all background regions in a given column (or
  row) will be fit with a polynomial of order "bkg_order" (see below).
  Values of "mean" and "median" compute the simple average and median,
  respectively, of the background region values in each column (or row).
  This parameter can also be specified in the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file. If
  specified in the reference file and given as an argument to the step by
  the user, the user-supplied value takes precedence.

``--bkg_order``
  The order of a polynomial function to be fit to the background
  regions.  The fit is done independently for each column (or row, if the
  dispersion is vertical) of the input image, and the fitted curve will be
  subtracted from the target data.  ``bkg_order`` = 0 (the minimum allowed
  value) means to fit a constant.  The user-supplied value (if any)
  overrides the value in the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file.  If neither is specified, a
  value of 0 will be used. If a sufficient number of valid data points is
  unavailable to construct the polynomial fit, the fit will be forced to
  0 for that particular column (or row). If "bkg_fit" is not "poly", this
  parameter will be ignored.

``--log_increment``
  For multi-integration extractions, if this parameter is set to a value greater
  than zero, an INFO-level log message will be printed every `log_increment` integrations
  to report on progress. Default value is 50.

``--save_profile``
  Flag to enable saving the spatial profile representing the extraction aperture or model.
  If True, the profile is saved to disk with suffix "profile".

``--save_scene_model``
  Flag to enable saving a model of the 2D flux as defined by the extraction aperture or PSF model.
  If True, the model is saved to disk with suffix "scene_model".

``--save_residual_image``
  Flag to enable saving the residual image (from the input minus the scene model)
  If True, the model is saved to disk with suffix "residual".

Step Arguments for IFU Data
---------------------------

``--center_xy``
  A list of two integer values giving the desired x/y location for the center
  of the circular extraction aperture used for extracting spectra from 3-D
  IFU cubes. Must be given in x,y order and in units of pixels along the x,y
  axes of the 3-D IFU cube, e.g. ``--center_xy="27,28"``.
  Default is None.

``--ifu_autocen``
  Switch to select whether or not to enable auto-centroiding of the extraction
  aperture for IFU point sources.  Auto-centroiding works by median collapsing the
  IFU cube across all wavelengths (shortward of 26 microns where the MRS throughput
  becomes extremely low) and using DAOStarFinder to locate the brightest
  source in the field. Default is ``False``.

``--bkg_sigma_clip``
  The background values will be sigma-clipped to remove outlier values from
  the determination of the background. The default value is a 3.0 sigma clip.

``--ifu_rfcorr``
  Switch to select whether or not to run 1d residual fringe correction on the
  extracted 1d spectrum (MIRI MRS only). Default is ``False``.

``--ifu_set_srctype``
  A string that can be used to override the extraction method for the source_type
  given by the SRCTYPE keyword. The allowed values are POINT and EXTENDED. The SRCTYPE keyword is
  not changed, instead the extraction method used is based on this parameter setting. This is
  only allowed for MIRI MRS IFU data. 

``--ifu_rscale``
   A float designating the number of PSF FWHMs to use for the extraction radius. This
   is a MIRI MRS only parameter. Values accepted are between 0.5 to 3.0. The default extraction
   size is set to 2 * FWHM. Values below 2 will result in a smaller
   radius, a value of 2 results in no change to radius and a value above 2 results in a larger
   extraction radius.

``--ifu_covar_scale``
   A float to be multiplied into the error arrays of the extracted spectra to account
   for covariance between adjacent spaxels in the IFU data cube.  The default value is
   1.0 (i.e., no correction) unless set by a user or a parameter reference file.  This
   parameter only affects MIRI and NIRSpec IFU spectroscopy.

Step Arguments for NIRISS SOSS Data
-----------------------------------

``--soss_atoca``
  Flag to enable using the ATOCA algorithm to treat order contamination. Default is ``True``.

``--soss_threshold``
  Threshold value for a pixel to be included when modeling the spectral trace. The default
  value is 0.01.

``--soss_n_os``
  An integer that sets
  the oversampling factor of the underlying wavelength grid used when modeling the
  trace. The default value is 2.

``--soss_wave_grid_in``
  Filename or SossWaveGridModel
  containing the wavelength grid used by ATOCA to model each valid pixel of the
  detector. If not given, the grid is determined based on an estimate of the flux
  (soss_estimate), the relative tolerance (soss_rtol) required on each pixel model
  and the maximum grid size (soss_max_grid_size).

``--soss_wave_grid_out``
  Filename to hold the wavelength
  grid calculated by ATOCA, stored in a SossWaveGridModel.

``--soss_estimate``
  Filename or SpecModel of the
  estimate of the target flux. The estimate must be a SpecModel with wavelength and
  flux values.

``--soss_rtol``
  The relative tolerance needed on a
  pixel model. It is used to determine the sampling of the soss_wave_grid when not
  directly given. Default value is 1.e-4.

``--soss_max_grid_size``
  The maximum grid size allowed. It is
  used when soss_wave_grid is not provided to make sure the computation time or the memory
  used stays reasonable. Default value is 20000.

``--soss_tikfac``
  This is the regularization
  factor used in the SOSS extraction. If not specified, ATOCA will calculate a
  best-fit value for the Tikhonov factor.

``--soss_width``
  This specifies the aperture
  width used to extract the 1D spectrum from the decontaminated trace. The default
  value is 40.0 pixels.

``--soss_bad_pix``
  This parameter sets the method
  used to handle bad pixels. There are currently two options: "model" will replace
  the bad pixel values with a modeled value, while "masking" will omit those pixels
  from the spectrum. The default value is "model".

``--soss_modelname``
  If set, this will provide
  the optional ATOCA model output of traces and pixel weights, with the filename
  set by this parameter. By default this is set to None and this output is
  not provided.
