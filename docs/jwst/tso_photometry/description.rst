Description
===========

:Class: `jwst.tso_photometry.TSOPhotometryStep`
:Alias: tso_photometry

The ``tso_photometry`` step does aperture photometry with a circular aperture
for the target.  Background is computed as the mean within a circular annulus.
The output is a table (ASCII ecsv format) containing the time at the
midpoint of each integration and the photometry values.

Assumptions
-----------
This step is intended to be used for aperture photometry with time-series
imaging exposures of a single source.  Only direct images should be used,
not spectra.

The input is usually assumed to be flux calibrated in units of MJy/sr, in which
case the output will be in Jy. If the flux calibration step was skipped, it
will convert from DN/s to the number of flat-fielded electrons measured per
integration.

The target is assumed to have been placed at the aperture reference location,
which is stored in the XREF_SCI and YREF_SCI FITS keywords
(note that these are 1-indexed values). Hence the step uses those keyword
values as the initial guess for the target location within the image.

Fit to Source
-------------

By default, before computing photometry, this step will attempt to refine the
planned position of the source by fitting a centroid to the data in a subimage
near the planned position.

For each integration, the source centroid is computed as the
center-of-mass for the subimage. The initial fit to the data uses a wide
search box centered on the planned position to derive an initial guess for
the centroid position. A secondary fit uses a smaller subimage to compute the
final centroid position. The subimage in both cases is background-subtracted
prior to the fit, from a median value of pixels outside an assumed source radius.
If the fit is successful, a Gaussian PSF is fit to the subimage at the
centroid location.

The centroid value, PSF 1-sigma width, and fit flux are reported in the catalog
for each integration. By default, however, the position used for computing
the aperture photometry is the median centroid, taken across the fit value for all
integrations.  If desired, the source position can be allowed to vary by
integration instead (set ``moving_centroid`` to True in the step arguments).

Aperture Photometry
-------------------
The Astropy affiliated package ``photutils`` does the work.

If the input file was *not* averaged over integrations (i.e. a _calints
product), and if the file contains an INT_TIMES table extension, the times
listed in the output table from this step will be extracted from column
'int_mid_MJD_UTC' of the INT_TIMES extension.  Otherwise,
the times will be computed from the exposure start time, the integration time,
and the number of integrations.  In either case, the times are
Modified Julian Date, time scale UTC.

If NaNs exist in the source or background annulus, they are masked and the value
returned is the sum over the real-valued data.
This is different from the convention for photometry in standard imaging mode:
in that case, a NaN value is returned if it is present in any of the pixels in
the source aperture. In the imaging case, the most likely cause of NaN pixels
in the aperture is NaN-valued central pixels for a saturated source in an image
with many sources. For TSO data, it is more likely that  single pixels are affected
in a given integration and science analysis will be focused on variability of one source.
If NaNs are present, the absolute flux will be underestimated, but the relative
values may still be useful.

Output Catalog
--------------

The output table contains these fields:

- ``MJD``: midpoint integration time
- ``aperture_sum``: source flux summed over the aperture
- ``aperture_sum_err``: error on the source flux
- ``annulus_sum``: background flux summed over the annular aperture
- ``annulus_sum_err`` : error on the background flux
- ``annulus_mean``: mean value in the background annulus (sum / area)
- ``annulus_mean_err``: error on the mean background value
- ``aperture_bkg``: background value (mean background * area)
- ``aperture_bkg_err``: error on the background value
- ``net_aperture_sum``: net source flux (summed flux minus background value)
- ``net_aperture_sum_err``: error on the net source flux
- ``aperture_x``: x center of the source aperture (zero-indexed)
- ``aperture_y``: y center of the source aperture (zero-indexed)
- ``aperture_ra_icrs``: RA position at the aperture center
- ``aperture_dec_icrs``: Dec position at the aperture center

If the source was centroided prior to aperture photometry (step argument
``centroid_source`` is True), the output table will additionally contain
these fields:

- ``centroid_x``: x centroid position fit to the source
- ``centroid_y``: y centroid position fit to the source
- ``psf_width_x``: x width (1-sigma) of a Gaussian PSF fit to the source
- ``psf_width_y``: y width (1-sigma) of a Gaussian PSF fit to the source
- ``psf_flux``: flux enclosed in a Gaussian PSF fit to the source


Subarrays
---------
If a subarray is used that is so small that the target aperture extends
beyond the limits of the subarray, the entire area of the subarray will be
used for the target aperture, and no background subtraction will be done.
A specific example is SUB64 with NIRCam, using PUPIL = WLP8.
