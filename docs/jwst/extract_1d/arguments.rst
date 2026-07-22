Step Arguments
==============
The ``extract_1d`` step has the following step-specific arguments.

The following arguments apply to **all modes** unless otherwise specified:

``--subtract_background`` (boolean, default=None)
  Flag to specify whether the background should be subtracted.  If None or `True`,
  background subtraction will be performed if there are background regions
  specified in the reference file.  If `False`, no background subtraction will be
  performed.  Has no effect for NIRISS SOSS data.

``--apply_apcorr`` (boolean, default=True)
  Switch to select whether or not to apply an APERTURE correction during
  processing. Has no effect for NIRISS SOSS data
  or for optimal extractions.

The following arguments apply to **slit and slitless spectroscopic data** only:

``--extraction_type`` (string, default="box")
  Specify the extraction type.

  * ``'box'``: A standard extraction is performed, summing over an aperture box.
    This is suitable for any input data (point sources and extended sources;
    resampled and unresampled images).
  * ``'optimal'``: A PSF-based optimal extraction is performed. This is
    best suited for un-resampled point sources. Currently, optimal extraction
    is only available for MIRI LRS Fixed Slit data.
  * None: Optimal extraction is attempted whenever ``use_source_posn`` is `True`.

``--use_source_posn`` (boolean, default=None)
  Specify whether the target and background extraction
  region locations specified in the :ref:`EXTRACT1D <extract1d_reffile>` reference
  file should be shifted to account for the expected position of the source. If None (the default),
  the step will decide whether to use the source position based
  on the observing mode and the source type. By default, source position corrections
  are attempted only for point sources in NIRSpec MOS/FS/BOTS and MIRI LRS fixed-slit exposures.
  Set to `False` to ignore position estimates for all modes; set to `True` to additionally attempt
  source position correction for extended sources.

``--position_offset`` (float, default=0)
  Specify a number of pixels (fractional pixels are allowed) to offset the
  extraction aperture from the nominal position.

``--model_nod_pair`` (boolean, default=True)
  Flag to enable fitting a negative trace during optimal extraction.
  If `True`, and ``extraction_type`` is ``'optimal'``, then a negative trace
  from nod subtraction is modeled alongside the positive source during
  extraction.  This will be attempted only if the input data has been background
  subtracted and the dither pattern type indicates that 2 nods were used.

``--optimize_psf_location`` (boolean, default=True)
  Flag to enable PSF location optimization during optimal extraction.
  If `True`, and ``extraction_type`` is ``'optimal'``, then the placement of
  the PSF model for the source location (and negative nod, if present)
  will be iteratively optimized. This parameter is recommended if
  negative nods are modeled.

``--smoothing_length`` (integer, default=None)
  If this value greater than 1 (and is an odd integer), the
  image data used to perform background extraction will be smoothed in the
  dispersion direction with a boxcar of this width.  If ``smoothing_length``
  is None (the default), the step will attempt to read the value from the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file.
  If a value is specified in the reference file,
  that value will be used.  Note that by specifying this parameter in the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file, a different value can
  be designated for each slit, if desired.
  If no value is specified either by the user nor in the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file,
  no background smoothing is done.

``--bkg_fit`` (string, default=None)
  The type of fit to perform to the background data in each image column
  (or row, if the dispersion is vertical). There are four allowed values:
  "poly", "mean", "median", and None (the default value). If None,
  the step will search the reference file for a value - if none is found,
  ``bkg_fit`` will be set to "poly". If set to "poly", the background
  values for each pixel within all background regions in a given column (or
  row) will be fit with a polynomial of order ``bkg_order`` (see below).
  Values of "mean" and "median" compute the simple average and median,
  respectively, of the background region values in each column (or row).
  This parameter can also be specified in the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file. If
  specified in the reference file and given as an argument to the step by
  the user, the user-supplied value takes precedence.

``--bkg_order`` (integer, default=None)
  The order of a polynomial function to be fit to the background
  regions.  The fit is done independently for each column (or row, if the
  dispersion is vertical) of the input image, and the fitted curve will be
  subtracted from the target data.  ``bkg_order=0`` (the minimum allowed
  value) means to fit a constant.  The user-supplied value (if any)
  overrides the value in the
  :ref:`EXTRACT1D <extract1d_reffile>` reference file.  If neither is specified, a
  value of 0 will be used. If a sufficient number of valid data points is
  unavailable to construct the polynomial fit, the fit will be forced to
  0 for that particular column (or row). If ``bkg_fit`` is not "poly", this
  parameter will be ignored.

``--log_increment`` (integer, default=50)
  For multi-integration extractions, if this parameter is set to a value greater
  than zero, an INFO-level log message will be printed every ``log_increment`` integrations
  to report on progress.

``--save_profile`` (boolean, default=False)
  Flag to enable saving the spatial profile representing the extraction aperture or model.
  If `True`, the profile is saved to disk with suffix "profile".

``--save_scene_model`` (boolean, default=False)
  Flag to enable saving a model of the 2D flux as defined by the extraction aperture or PSF model.
  If `True`, the model is saved to disk with suffix "scene_model".

``--save_residual_image``(boolean, default=False)
  Flag to enable saving the residual image (from the input minus the scene model).
  If `True`, the model is saved to disk with suffix "residual".

The following arguments apply to **IFU data** only:

``--center_xy`` (string, default=None)
  A string representing a list of two integer values giving the desired x/y location for the center
  of the circular extraction aperture used for extracting spectra from 3-D
  IFU cubes. Must be given in x,y order and in units of pixels along the x,y
  axes of the 3-D IFU cube, e.g., ``--center_xy="27,28"``.

``--ifu_autocen`` (boolean, default=False)
  Switch to select whether or not to enable auto-centroiding of the extraction
  aperture for IFU point sources.  Auto-centroiding works by median collapsing the
  IFU cube across all wavelengths (shortward of 26 microns where the MRS throughput
  becomes extremely low) and using `~photutils.detection.DAOStarFinder` to locate the brightest
  source in the field.

``--bkg_sigma_clip`` (float, default=3.0)
  The background values will be sigma-clipped to remove outlier values from
  the determination of the background.

``--ifu_rfcorr`` (boolean, default=True)
  Switch to select whether or not to run 1D residual fringe correction on the
  extracted 1D spectrum (MIRI MRS only).

``--ifu_set_srctype`` (string, default=None)
  A string that can be used to override the extraction method for the source type
  given by the SRCTYPE keyword. The allowed values are ``'POINT'`` and ``'EXTENDED'``. The SRCTYPE keyword is
  not changed, instead the extraction method used is based on this parameter setting. This is
  only allowed for MIRI MRS IFU data.

``--ifu_rscale`` (float, default=None)
   The number of PSF FWHMs to use for the extraction radius. This
   is a MIRI MRS only parameter. Values accepted are between 0.5 to 3.0. The default extraction
   size is set to ``2 * FWHM``. A value below 2 will result in a smaller
   radius, of 2 results in no change to radius, and above 2 results in a larger
   extraction radius.

``--ifu_covar_scale`` (float, default=1.0)
   A scale to be multiplied into the error arrays of the extracted spectra to account
   for covariance between adjacent spaxels in the IFU data cube.  The default value of
   1.0 means no correction. It can be set by a user or a parameter reference file.  This
   parameter only affects MIRI and NIRSpec IFU spectroscopy.

The following arguments apply to **NIRISS SOSS data** only:

``--soss_atoca`` (boolean, default=True)
  Flag to enable using the ATOCA algorithm to treat order contamination.

``--soss_threshold`` (float, default=1e-2)
  Threshold value for a pixel to be included when modeling the spectral trace.

``--soss_n_os`` (integer, default=2)
  The oversampling factor of the underlying wavelength grid used when modeling the
  trace.

``--soss_wave_grid_in`` (string, default=None)
  Filename or `~stdatamodels.jwst.datamodels.SossWaveGridModel`
  containing the wavelength grid used by ATOCA to model each valid pixel of the
  detector. If not given, the grid is determined based on an estimate of the flux
  (``soss_estimate``), the relative tolerance (``soss_rtol``) required on each pixel model
  and the maximum grid size (``soss_max_grid_size``).

``--soss_wave_grid_out`` (string, default=None)
  Filename to hold the wavelength
  grid calculated by ATOCA, stored in a `~stdatamodels.jwst.datamodels.SossWaveGridModel`.

``--soss_estimate`` (string, default=None)
  Filename or `~stdatamodels.jwst.datamodels.SpecModel` of the
  estimate of the target flux. The estimate must be a
  `~stdatamodels.jwst.datamodels.SpecModel` with wavelength and
  flux values.

``--soss_rtol`` (float, default=1e-4)
  The relative tolerance needed on a pixel model.
  It is used to determine the sampling of the SOSS wavelength grid when not
  directly given.

``--soss_max_grid_size`` (integer, default=20000)
  The maximum grid size allowed. It is used when SOSS wavelength grid is not
  provided to make sure the computation time or the memory
  used stays reasonable.

``--soss_tikfac`` (float, default=None)
  This is the regularization factor used in the SOSS extraction.
  If not specified, ATOCA will calculate a best-fit value for the Tikhonov factor.

``--soss_width`` (float, default=40)
  The aperture width (in pixels) used to extract the 1D spectrum
  from the decontaminated trace.

``--soss_bad_pix`` (string, default="masking")
  This parameter sets the method used to handle bad pixels.
  There are currently two options:

  * "model" will replace the bad pixel values with a modeled value, while
  * "masking" will omit those pixels from the spectrum.

``--soss_modelname`` (string, default=None)
  Filename for the optional ATOCA model output of traces and pixel weights.
  By default (None), this output is not provided.

``--soss_order_3`` (boolean, default=True)
  Flag to enable including spectral order 3 in the extraction for SOSS.
