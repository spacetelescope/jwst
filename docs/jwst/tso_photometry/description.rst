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
exposures.  Only direct images should be used, not spectra.

The target is assumed to have been placed at the aperture reference location,
which is stored in the XREF_SCI and YREF_SCI FITS keywords
(note that these are 1-indexed values). Hence the step uses those keyword
values as the target location within the image.

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
