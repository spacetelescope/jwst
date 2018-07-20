Description
===========
The background subtraction step in the calibration pipeline performs
image-from-image subtraction in order to accomplish subtraction of background
signal. The step takes as input one target exposure, to which the
subtraction will be applied, and a list of one or more background exposures.
Two different approaches to background image subtraction are used, depending
on the observing mode. Imaging and most spectroscopic modes use one method,
while a special method is used for Wide-Field Slitless Spectroscopy (WFSS).
as described at the end of this section. What follows here is the method used
for all non-WFSS modes.

Non-WFSS Modes
--------------

If more than one background exposure is provided, they will be averaged
together before being applied to the target exposure. Iterative sigma clipping
is applied during the averaging process, to reject sources or other outliers.
The clipping is accomplished using the astropy.stats.sigma_clip function (see
http://http://docs.astropy.org/en/stable/api/astropy.stats.sigma_clip.html).
The background step allows users to supply values for the sigma_clip
parameters `sigma` and `maxiters` (see below), in order to control the
clipping operation.

The average background image is produced as follows:

 - Clip the combined SCI arrays of all background exposures
 - Compute the mean of the unclipped SCI values
 - Sum in quadrature the ERR arrays of all background exposures, clipping the
   same input values as determined from the SCI arrays, and convert the result
   to an uncertainty in the mean
 - Combine the DQ arrays of all background exposures using a bitwise-OR
   operation

The average background exposure is then subtracted from the target exposure.
The subtraction consists of the following operations:

 - The SCI array of the average background is subtracted from the SCI
   array of the target exposure

 - The ERR array of the target exposure is currently unchanged, until full
   error propagation is implemented in the entire pipeline

 - The DQ arrays of the average background and the target exposure are
   combined using a bitwise-OR operation


If the target exposure is a simple ImageModel, the background image is
subtracted from it. If the target exposure is in the form of a CubeModel
(e.g. the result of a time series exposure), the background image
is subtracted from each plane of the CubeModel.

WFSS Mode
---------

For Wide-Field Slitless Spectroscopy expsoures (NIS_WFSS or NRC_WFSS),
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
'COMPLETE' in the output product.

Step Arguments
==============
The background step has two optional arguments, both of which are used only
when the step is applied to non-WFSS exposures and are passed as arguments
to the ``sigma_clip`` function:

* ``--sigma``: The number of standard deviations to use for the clipping limit.
               Defaults to 3.

* ``--maxiters``: The number of clipping iterations to perform, or ``None`` to
                  clip until convergence is achieved. Defaults to ``None``.

Reference File
==============
The only mode for which this step uses reference files is Wide-Field
Slitless Spectroscopy.  For WFSS, there are two reference files.
Reference type "WFSSBKG" is the background reference image, and reference
type "WAVELENGTHRANGE" contains information about the range of wavelengths
in the exposure.  The latter is used, together with a source catalog, to
create a mask showing the locations of source spectra in the image, and
hence, where the background regions are.

CRDS Selection Criteria
-----------------------
WFSSBKG reference files are selected on the basis of INSTRUME, DETECTOR,
EXP_TYPE, FILTER, and PUPIL values in the input science data set.

WAVELENGTHRANGE reference files are selected on the basis of INSTRUME,
EXP_TYPE, PUPIL (NIRCam only), and MODULE (NIRCam only).

WFSSBKG Reference File Format
-----------------------------
WFSSBKG reference files are FITS files with 3 IMAGE extensions and 1 BINTABLE
extension. The FITS primary data array is assumed to be empty.
The characteristics of the three IMAGE extensions are as follows:

=======  =====  ==============  =========
EXTNAME  NAXIS  Dimensions      Data type
=======  =====  ==============  =========
SCI        2    ncols x nrows   float
ERR        2    ncols x nrows   float
DQ         2    ncols x nrows   integer
=======  =====  ==============  =========

The BINTABLE extension contains the bit assignments used in the DQ array. It
uses ``EXTNAME=DQ_DEF`` and contains 4 columns:

* BIT: integer value giving the bit number, starting at zero
* VALUE: the equivalent base-10 integer value of BIT
* NAME: the string mnemonic name of the data quality condition
* DESCRIPTION: a string description of the condition

