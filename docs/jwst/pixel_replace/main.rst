Description
===========

:Classes: `jwst.pixel_replace.pixel_replace_step.PixelReplaceStep`
:Alias: pixel_replace

During 1-D spectral extraction (:ref:`extract_1d <extract_1d_step>` step),
pixels flagged as bad are ignored in the summation process.
If a bad pixel is part of the point-spread function (PSF) at a given wavelength, the
absence of the signal in the flagged pixel will lead to a hollow space at that wavelength in
the extracted spectrum.

For MIRI MRS and NIRSpec IFU exposures, NaN-valued pixels are partially compensated for during the
IFU cube building process by the overlap between detector pixels and output cube voxels.
The effects of NaN values are thus not as severe as for slit spectra, but can manifest
as small dips in the extracted spectrum when a NaN value lands atop the peak of a spectral
trace and cube building interpolates from lower-flux adjacent values.

To avoid these defects in the 1-D spectrum, this step estimates the flux values of pixels
flagged as ``DO_NOT_USE`` in 2-D extracted spectra using interpolation methods,
prior to rectification in the :ref:`resample_spec <resample_spec_step>` or
:ref:`cube_build <cube_build_step>` steps.
The ``pixel_replace`` step inserts these estimates into the 2-D spectral image array,
unsets the ``DO_NOT_USE`` flag in the DQ array and sets the ``FLUX_ESTIMATED`` flag for
each affected pixel instead. Error values and variance components for the replaced pixels
are similarly updated with estimated values, following the same interpolation method as is
used for the data.

This step is provided as a cosmetic feature and, for that reason, should be used with caution.

Algorithms
----------

Minimum Gradient Estimator
^^^^^^^^^^^^^^^^^^^^^^^^^^

The minimum gradient estimator (``algorithm = "mingrad"``) is the default algorithm.
It uses entirely local information to fill in missing pixel values.

This method tests the gradient along the spatial and spectral axes using immediately adjacent
pixels.  It chooses whichever dimension has the minimum absolute gradient and replaces the missing
pixel with the average of the two adjacent pixels along that dimension.  Near point sources,
this will thus favor replacement along the spectral axis due to spatial undersampling of the
PSF profile, while near bright extended emission lines it will favor replacement along the
spatial axis due to the steep spectral profile.

No replacement is attempted if a NaN value is bordered by another NaN value along a given axis.


Trace Modeling
^^^^^^^^^^^^^^

For some exposure types, it is possible to build a detailed spectral trace model
with the :ref:`adaptive_trace_model step <adaptive_trace_model_step>`.  The trace model
usually includes estimated flux values only for high signal-to-noise regions over the core
of the PSF, but it is generally a very good estimate for those values.
If ``algorithm = "trace_model"``, the ``pixel_replace`` step will preferentially replace
any missing values with estimates from the trace model if present.  Any remaining missing
values will be interpolated with the minimum gradient method, above.

If the :ref:`adaptive_trace_model step <adaptive_trace_model_step>` was not run prior to calling
``pixel_replace``, then it will be run at the start of the step, with the step parameter
``oversample=1`` and otherwise default values.

The errors and variances associated with pixels replaced from the trace model must still
be interpolated, since the trace model does not provide error estimates.  For these pixels,
the error values are linearly interpolated from the nearest valid pixels along the dispersion
direction.

For time series observations (TSO), the trace model is a 2-D image, generated from a median
spectral image across integrations.  Pixels replaced from the trace model for TSO observations
will have the same value in every integration.


Adjacent Profile Approximation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another option for interpolation is a method that uses a local spatial profile to interpolate
missing values (``algorithm = "fit_profile"``).

First, the input 2-D spectral cutout is scanned across the dispersion axis to determine
which cross-dispersion vectors (column or row, depending on dispersion direction) contain
at least one flagged pixel. Next, for each affected vector, a median normalized profile is created.

The adjacent arrays (the number of which is set by the step argument
``n_adjacent_cols``) are individually normalized. Next, each pixel in the profile is set to
the median of the normalized values. This results in a median of normalized values filling the vector.

Finally, this profile is scaled to the vector containing a missing pixel, and the value is
estimated from the scaled profile.

For some exposure types, PSF undersampling combined with the curvature of spectral traces on the
detector can lead the model-based adjacent profile estimator to derive incorrect values in the
vicinity of emission lines.  In these cases, either the minimum gradient or trace modeling methods
are likely to be more effective and reliable.

Reference Files
---------------
The ``pixel_replace`` step does not use any reference files.
