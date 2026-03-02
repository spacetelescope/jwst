Description
===========

:Class: `jwst.background.background_step.BackgroundStep`
:Alias: bkg_subtract

The background subtraction step performs
image-from-image subtraction in order to accomplish subtraction of background
signal.

Two different approaches to background image subtraction are used, depending
on the observing mode. Imaging and most spectroscopic modes use one method,
while a special method is used for Wide-Field Slitless Spectroscopy (WFSS).
For both background methods the output results are always returned in a new
data model, leaving the original input model unchanged.

Upon successful completion of the step, the ``S_BKDSUB`` keyword will be set to
"COMPLETE" in the output product.

This type of background subtraction is just one method available within the
JWST pipeline. See :ref:`Background Subtraction Methods <background_subtraction_methods>`
for an overview of all the methods and to which observing modes they are
applicable.

Input Data
----------
When run as part of a Stage 2 pipeline, the background step takes as input a
target exposure, to which the subtraction will be applied, and, for non-WFSS
modes, a list of one or more background exposures.  For WFSS modes, background
exposures are not required, but the target exposure must have a WCS applied,
and it must have a source catalog identified prior to processing.

When run as a standalone step, the input may be specified either as filenames
for appropriate intermediate products or as a filename for an association
that contains these filenames.

For example, for a non-WFSS mode observation, the background subtraction step can be run
on an intermediate target exposure ``exp_001``, subtracting ``exp_002`` and ``exp_003`` rate files
as follows::

    strun bkg_subtract exp_001_assignwcsstep.fits --bkg_list=exp_002_rate.fits,exp_003_rate.fits

For a WFSS mode observation, the step can be run on a single file, if it contains an assigned
WCS and the source catalog name in the ``SCATFILE`` FITS header keyword.  It can also be run on an
association file that contains the intermediate exposure filename and a source catalog.
For example, if the following contents are stored in a file called ``bkg_sub_asn.json``::

        {"products": [
            {
            "name": "wfss_exp_001",
            "members": [
                {
                    "expname": "wfss_exp_001_assignwcsstep.fits",
                    "exptype": "science",
                },
                {
                    "expname": "image_exp_001_cat.ecsv",
                    "exptype": "sourcecat"
                }
            ]
        }

then the background step can be run on the target exposure ``wfss_exp_001`` as::

    strun bkg_subtract bkg_sub_asn.json

From Python code, the source catalog can also be directly specified in the input model
before calling the BackgroundStep; e.g.::

    model = datamodels.open("wfss_exp_001_assignwcsstep.fits")
    model.meta.source_catalog = "image_exp_001_cat.ecsv"
    result = BackgroundStep.call(model)


Imaging and Non-WFSS, Non-SOSS Spectroscopic Modes
--------------------------------------------------
If more than one background exposure is provided, they will be averaged
together before being subtracted from the target exposure. Iterative sigma
clipping is applied during the averaging process, to reject sources or other
outliers.
The clipping is accomplished using :func:`astropy.stats.sigma_clip`.
The background step allows users to supply values for the :func:`astropy.stats.sigma_clip`
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
subtracted from it. If the target exposure is in the form of a 3-D
`~stdatamodels.jwst.datamodels.CubeModel`
(e.g., the result of a time series exposure), the average background image
is subtracted from each plane of the `~stdatamodels.jwst.datamodels.CubeModel`.

The combined, averaged background image can be saved using the step parameter
``save_combined_background``.

WFSS Mode
---------
For Wide-Field Slitless Spectroscopy exposures,
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
by the primary header keyword ``SCATFILE``), in conjunction with a reference file
that gives the wavelength range (based on filter and grism) that is relevant
to the WFSS image. All regions of the image that are free of source spectra
are used for scaling the background reference image.

A background mask is created and set to `True` where there are no sources, i.e., regions
where the background can be used.
This mask will be saved in the ``MASK`` extension of the intermediate output
file, saved with suffix "bsub", and will be accessible in the ``mask`` attribute of the
output datamodel.

The step argument ``wfss_mmag_extract`` can be used, if
desired, to set the minimum (faintest) AB mag of the source catalog objects
used to define the background regions.
The default is to use all source catalog entries that result in a spectrum
falling within the WFSS image.

The step argument ``wfss_mask`` can be used to provide a custom user mask
that overrides the source-catalog-derived mask. The argument should point to
a FITS or ASDF file openable as `~stdatamodels.jwst.datamodels.ImageModel`
containing a 2D array of integers in its ``.mask`` attribute (FITS ``MASK`` extension)
with pixels to be used as background set to 1 and other pixels set to 0.
The output of a previous run of the background subtraction step,
saved with suffix "bsub", can be used as such a custom mask
(after editing the ``.mask`` attribute as desired). To generate
a valid model from scratch, use something like::

    from stdatamodels.jwst.datamodels import ImageModel
    import numpy as np

    mask_data = np.ones((2048, 2048), dtype=np.uint32)  # same shape as science data
    mask_data[500:1500, 500:1500] = 0  # example: mask out central region
    mask_model = ImageModel()
    mask_model.mask = mask_data
    mask_model.save('custom_mask.fits')

To apply the user-defined mask as-is, the ``wfss_maxiter`` argument should be set to 0;
doing so will prevent any additional outlier rejection from being applied by the step.

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
