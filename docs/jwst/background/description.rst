Description
===========
The background subtraction step performs
image-from-image subtraction in order to accomplish subtraction of background
signal. The step takes as input one target exposure, to which the
subtraction will be applied, and a list of one or more background exposures.
Two different approaches to background image subtraction are used, depending
on the observing mode. Imaging and most spectroscopic modes use one method,
while a special method is used for Wide-Field Slitless Spectroscopy (WFSS).

Non-WFSS Modes
--------------
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

The average background image is produced as follows:

 * Clip the combined SCI arrays of all background exposures
 * Compute the mean of the unclipped SCI values
 * Sum in quadrature the ERR arrays of all background exposures, clipping the
   same input values as determined for the SCI arrays, and convert the result
   to an uncertainty in the mean
 * Combine the DQ arrays of all background exposures using a bitwise OR
   operation

The average background exposure is then subtracted from the target exposure.
The subtraction consists of the following operations:

 * The SCI array of the average background is subtracted from the SCI
   array of the target exposure

 * The ERR array of the target exposure is currently unchanged, until full
   error propagation is implemented in the entire pipeline

 * The DQ arrays of the average background and the target exposure are
   combined using a bitwise OR operation

If the target exposure is a simple ImageModel, the background image is
subtracted from it. If the target exposure is in the form of a 3-D CubeModel
(e.g. the result of a time series exposure), the background image
is subtracted from each plane of the CubeModel.

WFSS Mode
---------
For Wide-Field Slitless Spectroscopy expsoures (NIS_WFSS and NRC_WFSS),
a background reference image is subtracted from the target exposure.
Before being subtracted, the background reference image is scaled to match the
signal level of the target data within background (source-free) regions of the
image. 

The locations of source spectra are determined from a source catalog (specified
by the primary header keyword SCATFILE), in conjunction with a reference file
that gives the wavelength range (based on filter and grism) that is relavant
to the target data. All regions of the image that are free of source spectra
are used for scaling the background reference image. Robust mean values are
obtained for the background regions in the target image and for the same
regions in the background reference image, and the ratio of those two mean
values is used to scale the background reference image. The robust mean is
computed by excluding the lowest 25% and highest 25% of the data (using the
numpy.percentile function), and taking a simple arithmetic mean of the
remaining values.  Note that NaN values (if any) in the background
reference image are currently set to zero.  If there are a lot of NaNs,
it may be that more than 25% of the lowest values will need to be excluded.

For both background methods the output results are always returned in a new
data model, leaving the original input model unchanged.

Upon successful completion of the step, the S_BKDSUB keyword will be set to
"COMPLETE" in the output product.

