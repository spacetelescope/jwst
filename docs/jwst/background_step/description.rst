Description
===========

:Class: `jwst.background.BackgroundStep`
:Alias: bkg_subtract

The background subtraction step performs
image-from-image subtraction in order to accomplish subtraction of background
signal. The step takes as input a Spec2 association file or one target exposure,
to which the subtraction will be applied, and a list of one or more
background exposures.

Two different approaches to background image subtraction are used, depending
on the observing mode. Imaging and most spectroscopic modes use one method,
while a special method is used for Wide-Field Slitless Spectroscopy (WFSS).

This type of background subtraction is just one method available within the
JWST pipeline. See :ref:`Background Subtraction <background_subtraction>`
for an overview of all the methods and to which observing modes they're
applicable.

Imaging and Non-WFSS, Non-SOSS Spectroscopic Modes
--------------------------------------------------
If more than one background exposure is provided, they will be averaged
together before being subtracted from the target exposure. Iterative sigma
clipping is applied during the averaging process, to reject sources or other
outliers.
The clipping is accomplished using the function
`astropy.stats.sigma_clip
<http://docs.astropy.org/en/stable/api/astropy.stats.sigma_clip.html>`_.
The background step allows users to supply values for the ``sigma_clip``
parameters ``sigma`` and ``maxiters`` (see :ref:`bkg_step_args`),
in order to control the clipping operation.

For imaging mode observations, the calculation of the average background
image depends on whether the background exposures are "rate" (2D) or
"rateint" (3D) exposures. In the case of "rate" exposures, the average
background image is produced as follows:

#. Clip the combined SCI arrays of all background exposures. For mixtures
   of full chip and subarray data, only overlapping regions are used
#. Compute the mean of the unclipped SCI values
#. Sum in quadrature the ERR arrays of all background exposures, clipping the
   same input values as determined for the SCI arrays, and convert the result
   to an uncertainty in the mean
#. Combine the DQ arrays of all background exposures using a bitwise OR
   operation

In the case of "rateint" exposures, each background exposure can have multiple
integrations, so calculations are slightly more involved. The "overall" average
background image is produced as follows:

#. Clip the SCI arrays of each background exposure along its integrations
#. Compute the mean of the unclipped SCI values to yield an average image for
   each background exposure
#. Clip the means of all background exposure averages
#. Compute the mean of the unclipped background exposure averages to yield the
   "overall" average background image
#. Sum in quadrature the ERR arrays of all background exposures, clipping the
   same input values as determined for the SCI arrays, and convert the result
   to an uncertainty in the mean (This is not yet implemented)
#. Combine the DQ arrays of all background exposures, by first using a bitwise
   OR operation over all integrations in each exposure, followed by doing by a
   bitwise OR operation over all exposures.
        
The average background exposure is then subtracted from the target exposure.
The subtraction consists of the following operations:

#. The SCI array of the average background is subtracted from the SCI
   array of the target exposure

#. The ERR array of the target exposure is currently unchanged, until full
   error propagation is implemented in the entire pipeline

#. The DQ arrays of the average background and the target exposure are
   combined using a bitwise OR operation

If the target exposure is a simple ImageModel, the background image is
subtracted from it. If the target exposure is in the form of a 3-D CubeModel
(e.g. the result of a time series exposure), the average background image
is subtracted from each plane of the CubeModel.

The combined, averaged background image can be saved using the step parameter
``save_combined_background``.

WFSS Mode
---------
For Wide-Field Slitless Spectroscopy expsoures (NIS_WFSS and NRC_WFSS),
a background reference image is subtracted from the target exposure.
Before being subtracted, the background reference image is scaled to match the
signal level of the WFSS image within background (source-free) regions of the
image. The scaling factor is determined based on the variance-weighted mean
of the science data, i.e., ``factor = sum(sci*bkg/var) / sum(bkg*bkg/var)``.
This factor is equivalent to solving for the scaling constant applied to the
reference background that gives the maximum likelihood of matching 
the science data.
Outliers are rejected iteratively during determination of the scaling factor
in order to avoid biasing the scaling factor based on outliers. The iterative
rejection process is controlled by the
``wfss_outlier_percent``, ``wfss_rms_stop``, and ``wfss_maxiter`` step arguments.

The locations of source spectra are determined from a source catalog (specified
by the primary header keyword SCATFILE), in conjunction with a reference file
that gives the wavelength range (based on filter and grism) that is relevant
to the WFSS image. All regions of the image that are free of source spectra
are used for scaling the background reference image. The step argument
``wfss_mmag_extract`` can be used, if desired, to set the minimum (faintest)
abmag of the source catalog objects used to define the background regions.
The default is to use all source catalog entries that result in a spectrum
falling within the WFSS image.

For both background methods the output results are always returned in a new
data model, leaving the original input model unchanged.

Upon successful completion of the step, the S_BKDSUB keyword will be set to
"COMPLETE" in the output product.

SOSS Mode
---------
In a similar manner to WFSS modes, the NIRISS SOSS mode uses a set of reference
background templates, primarily for removal of flux contribution from zodiacal
dust.

First, a mask is derived to determine which regions of the input science data are
relatively uncontaminated, using a cutoff on flux percentile to mask out bright
regions of the integration. Then the mask is split into two components, one for
either side of a discontinuity in the SOSS background levels, a result of
instrumental effects. The mask on the right side of the detector is truncated
at column 950; pixels right of this column were found to lower the fitting accuracy
regardless of flux cutoff. The step then performs a best-fit analysis by scaling
each template in the background reference file to the data and finding the minimum
residual RMS error in the fitted background pixels. The best-fit template
is used to calculate and subtract the background for the entire science array.