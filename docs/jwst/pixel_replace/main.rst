Description
===========

:Classes: `jwst.pixel_replace.PixelReplaceStep`
:Alias: pixel_replace

During spectral extraction, pixels flagged as bad are ignored in the summation process.
If a bad pixel is part of the point-spread function (PSF) at a given wavelength, the
absence of the signal in the flagged pixel will lead to a divot at that wavelength in
the extracted spectrum.

To avoid this defect in the 1-D spectrum, this step estimates the flux values of pixels
flagged as ``DO_NOT_USE`` in 2-D extracted spectra, prior to rectification in the
``resample_spec`` step. ``pixel_replace`` inserts these estimates into the data array,
unsets the ``DO_NOT_USE`` flag and sets the ``FLUX_ESTIMATED`` flag for each pixel affected.

Algorithms
==========

Currently, one algorithm has been tested to estimate the missing fluxes.

Adjacent Profile Approximation
------------------------------

First, the input 2-d spectral cutout is scanned across the dispersion axis to determine
which cross-dispersion vectors (column or row, depending on dispersion direction) contain
at least one flagged pixel. Next, for each affected vector, a median normalized profile is created.

First, the adjacent arrays (the number of which is set by the step argument
``n_adjacent_cols``) are individually normalized. Next, each pixel in the profile is set to
the median of the normalized values. This results in a median of normalized values filling the vector.

Finally, this profile is scaled to the vector containing a missing pixel, and the value is
estimated from the scaled profile.