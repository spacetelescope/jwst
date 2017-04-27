Description
============

This step determines the mean count rate for each pixel by performing a linear
fit to the data in the input (jump) file.  The fit is done using "ordinary
least squares" (the "generalized least squares" is no longer an option).
The fit is performed independently for each pixel.  There are up to three
output files. The primary output file (RampFit, giving the slope at each
pixel) is always produced.  If the input exposure contains more than one
integration, the resulting slope images from each integration are stored as a
data cube in a second output data product.  A third, optional output product
is also available and is produced only when the step parameter 'save_opt' is
True (the default is False).  The output values will be in units of counts (DN)
or counts per second.  Following a description of the fitting algorithm, these
three type of output files are detailed below.


The count rate for each pixel is determined by a linear fit to the
cosmic-ray-free ramp intervals for each pixel. CR-free intervals are derived
using the 4-D GROUPDQ array of the input data set, under the assumption that
the jump step will have already flagged CR's. Ramp intervals are also terminated
where saturation flags are found.  Ramp intervals that are noiseless, or have
no signal, or contain only 2 reads will initially have fits with variance = 0,
preventing their slopes from contributing to the weighted slopes.  In these
cases, the variance will be recalculated as the poisson noise of the ramp added
in quadrature to the pixel-specific read noise, ensuring that all variance
values are positive.  If the input dataset has only a single group in each
integration, the count rate for all unsaturated pixels in that integration will
be calculated to be the value of the science data in that group divided by the
exposure time.  If the input dataset has only two groups per integration, the
count rate for all unsaturated pixels in each integration will be calculated
from the 2 valid values of the science data.  If any input dataset contains
ramps saturated in their second read, the count rates for those pixels in that
integration will be calculated to be the value of the science data in that group
divided by the exposure time. After computing the slopes for all intervals for
a given pixel, the final slope is determined as a weighted average from all
intervals and is written to a file as the primary output product.  In this
output product, the 4-D GROUPDQ from all integrations is compressed into 2-D,
which is then merged (using a bitwise OR) with the input 2-D PIXELDQ to create
the output PIXELDQ.


If the input exposure contains more than one integration, the resulting slope
images from each integration are stored as a data cube in a second output data
product. Each plane of the 3-D SCI, ERR, and DQ arrays in this product is the
result for a given integration.  In this output product, the 4-D GROUPDQ from
an integration is compressed into 2-D, which is then merged with the input 2-D
PIXELDQ to create the output PIXELDQ for each integration.  By default, the
name of this output product is based on the name of the input file; the user
can override this name by specifying a name using the parameter int_name.


A third, optional output product is also available and is produced only when
the step parameter 'save_opt' is True (the default is False). This optional
product contains 4-D arrays called SLOPE, SIGSLOPE, YINT, SIGYINT, and WEIGHTS,
which contain the slopes, uncertainties in the slopes, y-intercept, uncertainty
in the y-intercept, and fitting weights for each ramp interval of each pixel.
The y-intercept refers to the result of the fit at an exposure time of zero.
This product also contains a 3-D array called PEDESTAL, which gives the signal
at zero exposure time for each pixel, and the 4-D CRMAG array, which contains
the magnitude of each read that was flagged as having a CR hit.  By default,
the name of this output file is based on the name of the input file; the user
can override this name by specifying a name using the parameter opt_name.  In
this optional output product, the pedestal array is calculated for each
integration by extrapolating the final slope (the weighted average of the slopes
of all of ramp segments in the integration) for each pixel from its value at
the first sample to an exposure time of zero. Any pixel that is saturated on
the first read is given a pedestal value of 0.  Before compression, the cosmic
ray magnitude array is equivalent to the input SCI array but with the only
nonzero values being those whose pixel locations are flagged in the input
GROUPDQ as cosmic ray hits. The array is compressed, removing all reads in which
all the values are 0 for pixels having at least one read with a non-zero
magnitude. The order of the cosmic rays within the ramp is preserved.


The fitting algorithm does an 'optimal' linear fit, which is the weighting used
by Fixsen et al, PASP,112, 1350. ('unweighted' in which pixels are equally
weighted, is no longer a weighting option.)  Pixels are processed simultaneously
in blocks using the array-based functionality of numpy.  The size of the block
depends on the image size and the number of groups.


Upon successful completion of this step, the status keyword S_RAMP will be set
to COMPLETE.


After the upcoming DMS build, the MIRI last frame correction step will be
updated to flag all the pixels in the last group of data in each integration
of a MIRI exposure, so that those data do not get used in either the jump
detection or ramp fitting steps.  As a result, the ramp fitting step will be
updated to make sure it handles the flagged group properly and does not include
any data from the last group of each integration in its calculations; for MIRI
exposures that have original values of 2 and 3 groups per integration, ramp
fitting processing will proceed using only the first 1 and 2 groups,
respectively, using the calculations described above.  For MIRI exposures that
have an original value of only 1 group per integration, the last group will NOT
be flagged by the last frame correction step, so that there will always be at
least 1 group of data to work with in subsequent steps.  Hence the special ramp
fitting processing that's applied to exposures that only have a single group
will be applied to MIRI exposures that originally have 1 and 2 groups.

Step Arguments
==============
The ramp fitting step has three optional arguments that can be set by the user:

* ``--save_opt``: A True/False value that specifies whether to write
  optional output information.

* ``--opt_name``: A string that can be used to override the default name
  for the optional output information.

* ``--int_name``: A string that can be used to override the default name
  for the integration-by-integration slopes, for the case that the input
  file contains more than one integration.


