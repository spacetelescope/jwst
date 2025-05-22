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
imaging exposures.  Only direct images should be used, not spectra.

The target is assumed to have been placed at the aperture reference location,
which is stored in the XREF_SCI and YREF_SCI FITS keywords
(note that these are 1-indexed values). Hence the step uses those keyword
values as the target location within the image.

The input is usually assumed to be flux calibrated in units of MJy/sr, in which
case the output will be in Jy. If the flux calibration step was skipped, it
will convert from DN/s to the number of flatfielded electrons measured per
integration.

Algorithm
---------
The Astropy affiliated package ``photutils`` does the work.

If the input file was *not* averaged over integrations (i.e. a _calints
product), and if the file contains an INT_TIMES table extension, the times
listed in the output table from this step will be extracted from column
'int_mid_MJD_UTC' of the INT_TIMES extension.  Otherwise,
the times will be computed from the exposure start time, the integration time,
and the number of integrations.  In either case, the times are
Modified Julian Date, time scale UTC.

If NaNs exist in the source or background annulus, they are masked out and the value
returned is the sum over the real-valued data.
This is different from the convention for photometry in standard imaging mode:
in that case, a NaN value is returned if it is present in any of the pixels in
the source aperture. In the imaging case, the mostly likely cause of NaN pixels
in the aperture is NaN-valued central pixels for a saturated source in an image
with many sources. For TSO data, it is more likely that  single pixels are affected
in a given integration and science analysis will be focused on variability of one source.
If NaNs are present, the absolute flux will be underestimated, but the relative values may still
be useful.


The output table contains these fields:

- MJD
- aperture_sum
- aperture_sum_err
- annulus_sum
- annulus_sum_err
- annulus_mean
- annulus_mean_err
- aperture_bkg
- aperture_bkg_err
- net_aperture_sum
- net_aperture_sum_err

Subarrays
---------
If a subarray is used that is so small that the target aperture extends
beyond the limits of the subarray, the entire area of the subarray will be
used for the target aperture, and no background subtraction will be done.
A specific example is SUB64 with NIRCam, using PUPIL = WLP8.
