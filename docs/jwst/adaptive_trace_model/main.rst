Description
===========

:Classes: `jwst.adaptive_trace_model.adaptive_trace_model_step.AdaptiveTraceModelStep`
:Alias: adaptive_trace_model


Overview
--------
The ``adaptive_trace_model`` step models the spectral trace in a 2D spectral image
with a set of univariate basis spline fits to the spatial profile, along the dispersion axis.
Optionally, the step may also use the trace model to oversample the input data by a
specified factor.

This step is intended in part to address
`spatial undersampling effects <https://jwst-docs.stsci.edu/known-issues/resampling-artifacts>`__
in NIRSpec IFU and MIRI MRS spectra extracted from rectified cubes.
This "resampling noise" manifests as low-frequency oscillations in spectra extracted from
apertures smaller than the observational point-spread-function (PSF). Interpolating the data
onto a higher resolution grid prior to building a rectified spectral cube can mitigate
these spectral artifacts.

NIRSpec fixed slit (FS) and MOS spectra and MIRI LRS slit spectra may also show similar
resampling artifacts due to undersampling that can be mitigated by oversampling prior
to creating the rectified 2D spectral image in the :ref:`resample_spec <resample_spec_step>`
step.

This step is currently available for all NIRSpec and MIRI spectral modes.
Trace models may be generated for time-series spectra, but oversampling is disabled
for this mode. The step is incorporated into the :ref:`calwebb_spec2 <calwebb_spec2>` and
:ref:`calwebb_spec3 <calwebb_spec3>` pipelines, prior to the
:ref:`pixel_replace <pixel_replace_step>` and :ref:`cube_build <cube_build_step>` or
:ref:`resample_spec <resample_spec_step>` steps. It may also be run as a standalone step.

Upon successful completion of this step, the status keyword S_TRCMDL is set to "COMPLETE",
and the trace model image is stored in the output datamodel's ``trace_model`` attribute
(FITS extension TRACEMODEL).  If oversampling was performed, the input flux, error, and variance
images are replaced with interpolated data sampled onto a new pixel grid. The data quality (DQ),
wavelength, and regions images are also oversampled, but any additional images
(e.g., pathloss corrections) are not propagated to the output datamodel.

Algorithm
---------
The adaptive spline modeling used by this step depends on two assumptions:
one, that the PSF should change only slowly with wavelength, and two, that unresolved point
sources should have centroids that remain fixed in celestial coordinates at all wavelengths.
These assumptions enable the creation of a model of the spatial profile at each wavelength
by fitting the spatial data within a window of each dispersion element. The data
are fit with a cubic basis spline function, with sufficient knots defined to evenly
sample the spatial profile without overfitting artifacts.

Model the Trace
^^^^^^^^^^^^^^^

Modeling is performed one spectral region (IFU slice or slit cutout) at a time.  Each
region is separately modeled as follows:

1. Determine if the region has sufficient signal to justify fitting spline models. For
   IFU, if the mean value for the slice is less than 10-sigma higher than the overall
   mean for the image, by default, then the slice is ignored and no modeling is performed.
   For slit-like modes, if the median signal-to-noise collapsed across wavelengths
   is greater than 10 (by default), the slit is ignored and no modeling is attempted.

#. Compute spatial (cross-dispersion) coordinates for every pixel in the region.

#. Make a normalized image by dividing by the sum over each dispersion element (column) in
   the region.  If significant negative artifacts are present, from subtracting nod pairs
   in the :ref:`bkg_subtract <background_subtraction>` step, then only positive data are
   summed.

#. For each column in the region:

   a. Select normalized data within a range of nearby wavelengths.

   #. Fit the normalized data by spatial coordinate with a cubic basis spline function.

   #. Reject any data points more than 2.5 sigma from the spline model and re-fit
      the remaining data, for a maximum of 3 iterations, by default. The threshold value and
      number of iterations may vary slightly by spectral mode.

   #. Evaluate the model at the input coordinates for the column and determine a scaling
      factor to reproduce the observed values, from the weighted mean ratio of the original
      fluxes to the normalized spline model.

The set of spline models and scale factors for each wavelength in each region constitutes
the adaptive trace model for the spectral image.

The spline modeling assumptions are generally only appropriate for compact sources, so an
additional check is made to determine regions for which the model is likely to be accurate.
The slope of the model flux is computed for each column pixel as the
absolute difference between the normalized spline model at that pixel and its immediate neighbor.
Slope values higher than a threshold value (step parameter ``slope_limit``) indicate
a compact source region.  The trace model will be evaluated for these regions, with some
padding for nearby pixels; it will not be evaluated in other regions.

For MOS and fixed slit data, note that the model accuracy also depends on the position of the
source within the spectral region.  Sources on the edge of a MOS slitlet, for example,
may not be well modeled by the spline fits, even if they are bright and compact.

If no oversampling is desired (i.e., the ``oversample`` parameter is set to 1.0), then the trace
model is evaluated at every input pixel in a compact source region to create a wavelength-dependent
spatial profile.  This image is stored in the output datamodel, in the ``trace_model`` attribute.
Regions for which a spline model could not be computed, or which did not meet the compact source
criteria, are set to NaN in the image. The step then returns without further changes to the input
datamodel.  The rest of the algorithm description, below, applies only to oversampling.

For IFU data, the trace modeling portion of the code can optionally use multiprocessing
to improve runtime by modeling the slices independently in separate processes.

Oversample the Flux
^^^^^^^^^^^^^^^^^^^

If oversampling is desired, the step will create new data arrays. The spatial dimension is scaled
by the oversampling factor; the spectral dimension remains the same.  Each region is again
processed separately, by interpolating the spectral flux onto the new grid as follows:

1. Compute spatial coordinates for each pixel in the oversampled grid.

#. For each column in the region:

   a. Compute a linear interpolation of the data onto the new oversampled coordinates (:math:`f_{linear}`).

   #. If a spline fit is available for this column, evaluate the spline model at the original
      coordinates for the column.

   #. Construct the residual between the evaluated spline fit and the original data, and
      linearly interpolate it onto the oversampled coordinates (:math:`f_{residual}`).

   #. Compute the slope of each column pixel as the absolute difference between the normalized
      spline model at that pixel and its immediate neighbor.

   #. Evaluate the spline model at the oversampled coordinates (:math:`f_{spline}`).

#. Construct the oversampled region flux (*f*) from a piecewise model:

      :math:`f = f_{spline} + f_{residual}`, where the slope is greater than a threshold value

   and

      :math:`f = f_{linear}`, otherwise.

This method results in an interpolated flux image that uses the spline models for any bright,
compact sources and a linear interpolation for faint, diffuse regions. The residual image added
into the spline model accounts for any local structures that are not well modeled by the spline
profile.

Note that the interpolation process may provide output values for some pixels corresponding to
data with NaN values in the input.  If the region is modeled by a valid spline interpolation, the
missing values are extrapolated and replaced with real values from the spline model plus residual
flux.  These values will be marked in the DQ plane with a FLUX_ESTIMATED flag (see below).

Optionally, if the ``psf_optimal`` step parameter is set to True, fit threshold and slope limits
are ignored, so that spline models are created and used for all pixels, and the residual image
is not added into the oversampled flux.  This option is only appropriate for simple, isolated
point sources, but if used, can significantly improve the signal-to-noise ratio (SNR) for
extracted spectra, at the cost of ignoring non-PSF structures.

Alternately, crowded fields with multiple stars may benefit particularly from setting
the ``fit_threshold`` and ``slope_limit`` parameters to zero in order to ensure proper
modeling of both bright and faint stars.

Alongside the oversampled flux image, the set of spline models evaluated at all compact source
coordinates (:math:`f_{spline}`, above, where the slope condition is met) are saved to the output
model in the ``trace_model`` attribute, as a record of the wavelength-dependent spatial profile
for the oversampled data.

Propagate DQ, Error, and Variance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To match the oversampled flux image, the error and variance arrays in the datamodel are
linearly interpolated onto the oversampled grid.  The DQ array is oversampled as well,
with a nearest-pixel interpolation.  It is then updated with a FLUX_ESTIMATED flag for any
pixels that were NaN in the input but were replaced with real values by the spline modeling.

After oversampling, error and variance arrays are inflated by a factor dependent on the
oversampling ratio to account for the introduced covariance.  This factor is intended
to produce resampled cubes or images with the same SNR for all oversampling factors, to first order.
The value of the factor (*X*) was empirically determined from tests on a line-free region
of a stellar spectrum, and is calculated for oversample factor *N* as:

.. math::
   X = 0.23 N + 0.77

The oversampled error image is multiplied by X; variance images are multiplied
by X\ :sup:`2`.

Note that the inflated error arrays do not accurately reflect the per-pixel errors
on the oversampled flux, but rather are intended to produce approximately correct
errors after further resampling. The oversampled product is primarily intended to
be an intermediate format, prior to building a rectified spectral cube or image.

Update the WCS
^^^^^^^^^^^^^^

Finally, for oversampled data, the WCS object for the exposure must be updated to include a
transform from the oversampled pixel coordinates to the original detector coordinates.  The transform
is stored in a frame called "coordinates", prior to the "detector" frame.  WCS operations following
oversampling should use "coordinates" as the input frame.  For example, to retrieve the
world coordinates for pixel x, y in the oversampled image, these operations are equivalent::

   ra, dec, lam = oversampled_model.meta.wcs(x, y)

and::

   oversampled_transform = oversampled_model.meta.wcs.get_transform('coordinates', 'world')
   ra, dec, lam = oversampled_transform(x, y)

To retrieve the transform from original detector pixels to world coordinates instead,
use the "detector" frame::

   detector_pixel_transform = oversampled_model.meta.wcs.get_transform('detector', 'world')

Reference Files
---------------
The ``adaptive_trace_model`` step does not use any reference files.

References
----------

The adaptive trace model algorithm is based on work by D. Law and M. Clarke,
"Mitigating Resampling Artifacts for the JWST IFU Spectrometers with Adaptive Trace Modeling"
(`2026, AJ, 171, 304 <https://ui.adsabs.harvard.edu/abs/2026AJ....171..304L/abstract>`__).
