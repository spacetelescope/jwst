Description
===========

:Classes: `jwst.pixel_replace.pixel_replace_step.PixelReplaceStep`
:Alias: pixel_replace

During 1-D spectral extraction (:ref:`extract_1d <extract_1d_step>` step),
pixels flagged as bad are ignored in the summation process.
If a bad pixel is part of the point-spread function (PSF) at a given wavelength, the
absence of the signal in the flagged pixel will lead to a gap at that wavelength in
the extracted spectrum.

For MIRI MRS and NIRSpec IFU exposures, NaN-valued pixels are partially compensated for during the
IFU cube building process by the overlap between dithered detector pixels and output cube voxels.
The effects of NaN values are thus not as severe as for slit spectra, but can still be seen
when a NaN value lands atop the peak of a spectral trace. In that case, cube building reconstructs
the output voxel from lower-flux adjacent values, leading to a small dip in the extracted spectrum.

To minimize such artifacts in the 1-D spectra, this step uses interpolation methods to estimate the flux values of
missing pixels flagged as ``DO_NOT_USE`` in 2-D calibrated spectra
prior to rectification in the :ref:`resample_spec <resample_spec_step>` or
:ref:`cube_build <cube_build_step>` steps.
The ``pixel_replace`` step inserts these estimates into the 2-D spectral image array,
unsets the ``DO_NOT_USE`` flag in the DQ array, and sets the ``FLUX_ESTIMATED`` flag for
each affected pixel instead. Error values and variance components for the replaced pixels
are similarly updated with estimated values, following the same interpolation method as is
used for the data.

Empirically, this step can significantly reduce the number of artifacts in 1-D spectra by using information in
detector space to estimate the missing values.  However, any such interpolation should be treated with caution,
particularly where the input spectra are not smooth (i.e., in the vicinity of sharp spectral features).
Three different algorithms are provided for this step, each of which have their own benefits and limitations.


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

No replacement is attempted if a NaN value is bordered by another NaN value along a given axis (i.e. for anything
other than single isolated missing pixels).

This is the fastest method (typically seconds vs minutes per exposure), the least aggressive,
and generally performs well at interpolating
values for single missing pixels in most circumstances (point sources, extended sources, and emission
features).  However, it cannot estimate values for larger
blocks of bad pixels (such as 3x3 bad pixel regions often seen in NIRSpec data) and hence cannot fix some larger
spectral artifacts that the other two methods can.


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

The errors and variances associated with pixels replaced from the trace model must be
separately interpolated, since the trace model does not provide error estimates.  These error
values are linearly interpolated from the nearest valid pixels along the dispersion
direction.

For time series observations (TSO), the trace model is a 2-D image, generated from a median
spectral image across integrations.  Pixels replaced from the trace model for TSO observations
will have the same value in every integration.

This method can interpolate values for larger blocks of missing pixels near bright point sources using
a dynamically-defined model of the dispersed trace profile while retaining the performance of the "mingrad"
algorithm far from the center of the spectral trace.  However, it cannot interpolate values for larger blocks
of bad pixels in the wings of the PSF or in large extended sources.


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

This is the most aggressive pixel replacement method, which can interpolate values for larger blocks of missing
pixels across the entire detector, for both point and extended sources.
However, PSF undersampling combined with the curvature of spectral traces on the
detector can sometimes lead this estimator to derive incorrect values in the
vicinity of emission lines (for both isolated and large blocks of bad pixels).
In such cases, either the minimum gradient or trace modeling methods
are likely to be more effective and reliable.

Reference Files
---------------
The ``pixel_replace`` step does not use any reference files.
