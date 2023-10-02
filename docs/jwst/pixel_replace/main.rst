Description
===========

:Classes: `jwst.pixel_replace.PixelReplaceStep`
:Alias: pixel_replace

During 1-D spectral extraction (:ref:`extract_1d <extract_1d_step>` step),
pixels flagged as bad are ignored in the summation process.
If a bad pixel is part of the point-spread function (PSF) at a given wavelength, the
absence of the signal in the flagged pixel will lead to a hollow space at that wavelength in
the extracted spectrum.

To avoid this defect in the 1-D spectrum, this step estimates the flux values of pixels
flagged as ``DO_NOT_USE`` in 2-D extracted spectra using interpolation methods,
prior to rectification in the :ref:`resample_spec <resample_step>` step.
``pixel_replace`` inserts these estimates into the 2-D data array,
unsets the ``DO_NOT_USE`` flag, and sets the ``FLUX_ESTIMATED`` flag for each affected pixel.

This step is provided as a cosmetic feature and, for that reason, should be used with caution.

Algorithms
==========

Adjacent Profile Approximation
------------------------------

This is the default (and most extensively tested) algorithm for most spectroscopic modes.

First, the input 2-D spectral cutout is scanned across the dispersion axis to determine
which cross-dispersion vectors (column or row, depending on dispersion direction) contain
at least one flagged pixel. Next, for each affected vector, a median normalized profile is created.

The adjacent arrays (the number of which is set by the step argument
``n_adjacent_cols``) are individually normalized. Next, each pixel in the profile is set to
the median of the normalized values. This results in a median of normalized values filling the vector.

Finally, this profile is scaled to the vector containing a missing pixel, and the value is
estimated from the scaled profile.

Minimum Gradient Estimator
--------------------------

In the case of the MIRI MRS, NaN-valued pixels are partially compensated during the IFU cube building process
using the overlap between detector pixels and output cube voxels.  The effects of NaN values are thus not
as severe as for slit spectra, but can manifest as small dips in the extracted spectrum when a NaN value
lands atop the peak of a spectral trace and cube building interpolates from lower-flux adjacent values.

Pixel replacement can thus be useful in some science cases for the MIRI MRS as well, but undersampling combined with
the curvature of spectral traces on the detector
can lead the model-based adjacent profile estimator to derive incorrect values in the vicinity of
emission lines.  The minimum gradient estimator is thus another optional algorithm that uses entirely
local information to fill in the missing pixel values.

This method tests the gradient along the spatial and spectral axes using immediately adjacent pixels.  It chooses
whichever dimension has the minimum absolute gradient and replaces the missing pixel with the average of the
two adjacent pixels along that dimension.  Near point sources this will thus favor replacement along the spectral
axis due to spatial undersampling of the PSF profile, while near bright extended emission lines it will favor
replacement along the spatial axis due to the steep spectral profile.  No replacement is attempted if a NaN
value is bordered by another NaN value along a given axis.
