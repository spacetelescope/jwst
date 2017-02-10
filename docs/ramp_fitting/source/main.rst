Description
============

This step determines the mean count rate for each pixel by performing a linear
fit to the data in the input (jump) file.  The fit is performed independently
for each pixel.  The fit may be done using "ordinary least squares" or
"generalized least squares."  The primary output file (RampFit, giving the
slope at each pixel) is the same for both algorithms, and the per-integration
output file is also comparable, but the optional output file is significantly
different for generalized least squares.  The output values will be in
units of counts (DN) or counts per second.  Ordinary least squares is the
default and will be described first.

When using ordinary least squares, the count rate for each pixel is determined 
by a linear fit to the cosmic-ray-free ramp intervals for each pixel. CR-free 
intervals are derived using the 4-D GROUPDQ array of the input data set, under
the assumption that the jump step will have already flagged CR's. Ramp intervals
are also terminated where saturation flags are found.  Ramp intervals that are 
noiseless, or have no signal, or contain only 2 reads will initially have fits 
with variance = 0, preventing their slopes from contributing to the weighted 
slopes.  In these cases, the variance will be recalculated as the poisson noise 
of the ramp added in quadrature to the pixel-specific read noise, ensuring that 
all variance values are positive.  If the input dataset has only a single group
in each integration, the count rate for all unsaturated pixels in that 
integration will be calculated to be the value of the science data in that group 
divided by the exposure time.  If the input dataset has only two groups per 
integration, the count rate for all unsaturated pixels in each integration 
will be calculated from the 2 valid values of the science data.  If any input 
dataset contains ramps saturated in their second read, the count rates for those 
pixels in that integration will be calculated to be the value of the science 
data in that group divided by the exposure time. After computing the slopes 
for all intervals for a given pixel, the final slope is determined as a weighted 
average from all intervals and is written to a file as the primary output 
product.  In this output product, the 4-D GROUPDQ from all integrations is 
compressed into 2-D, which is then merged (using a bitwise OR) with the input 
2-D PIXELDQ to create the output PIXELDQ.  In addition, the output PIXELDQ array 
is updated to have DQ flags set to NO_GAIN_VALUE and DO_NOT_USE for all pixels 
that are non-positive or NaN in the gain array. 

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
This product also contains a 3-D array called PEDESTAL, which gives the
y-intercept at zero exposure time for each pixel, and the 4-D CRMAG array,
which contains the magnitude of each read that was flagged as having a CR hit.
By default, the name of this output file is based on the name of the input file;
the user can override this name by specifying a name using the parameter
opt_name.

In this optional output product, the pedestal array is calculated for each
integration by extrapolating the final slope for each pixel from its value at
the first sample to an exposure time of zero. Any pixel that is saturated on
the first read is given a pedestal value of 0.  Before compression, the cosmic
ray magnitude array is equivalent to the input SCI array but with the only
nonzero values being those whose pixel locations are flagged in the input
GROUPDQ as cosmic ray hits. The array is compressed, removing all reads in which
all the values are 0 for pixels having at least one read with a non-zero
magnitude. The order of the cosmic rays within the ramp is preserved.

By default, the OLS fitting algorithm does an unweighted (uniform) linear fit, 
meaning that all data points are equally weighted. Instead, the user can 
specify the parameter 'weighting' to be 'optimal', in which case fitting will 
use the weighting used by Fixsen et al, PASP,112, 1350.

The algorithm parameter can be used to select generalized least squares,
either by specifying the strun command-line option '--algorithm GLS', or by
including ``algorithm = "GLS"`` in the configuration file ramp_fit.cfg.  If
generalized least squares was specified, the fit will be done (for each
pixel) to the entire ramp, rather than to the CR-free segments, and the
amplitude of each jump as well as the overall slope will be computed.  It
is assumed that previous steps (jump and saturation) will have flagged
pixels in the GROUPDQ extension to indicate the pixels and reads that were
affected by cosmic-ray hits or saturation.  Since different pixels may have
been hit by different numbers of cosmic rays, pixels will be selected and
fit differently, depending on how many cosmic rays hits were detected in
the previous jump step.  For each pixel with no CR hits, a simple straight
line (intercept and slope) will be fit.  For each pixel with one CR hit,
the parameters of the fit will be intercept, slope, and the amplitude of a
jump at a known read (i.e., known from a flag in the GROUPDQ extension),
and similarly for pixels with multiple CR hits.  The fits will be done in
sets:  all pixels with no CR hits, all pixels with one CR hit, all pixels
with two CR hits, etc.  Within each such set, the function to be fit is
the same, except for the read at which the jumps occurred.  For each such
set, the problem is set up as a matrix equation of the form::

    Y = X * P

where Y is a (column) vector of the observed values of the ramp, in
electrons, i.e. each element corresponds to a read (or an average of a
group of reads, for grouped data); P is a vector of the parameters to be
determined (intercept, slope, and zero or more jump amplitudes); and X is
a matrix of (known) independent variables.  X has 2 + num_cr columns and
ngroups rows, where num_cr is the number of CR hits for pixels in the
current set, and ngroups is the number of groups in the ramp.  The first
column of X is all ones (for determining the intercept).  The second column
of X is the time at the end of each group (for finding the slope).  Each
subsequent column is a Heaviside function, zero for each row (i.e. for
each group) prior to the CR hit and one for each row on or after the CR
hit.  These columns are in time order of the CR hits; that is, the column
for the first CR hit would have fewer zeros and more ones than any
subsequent column.  For example, consider a ramp with a slope of one
electron per group and these values (including two CR hits, each with an
amplitude of 100 electrons):

====   =====
time   value
====   =====
10.7   1
21.4   2
32.1   3
42.8   104
53.5   105
64.2   206
74.9   207
85.6   208
====   =====

If the group time was 10.7 seconds, the X matrix would be:

===  ====   ==    ==
--------------------
1    10.7    0     0
1    21.4    0     0
1    32.1    0     0
1    42.8    1     0
1    53.5    1     0
1    64.2    1     1
1    74.9    1     1
1    85.6    1     1
===  ====   ==    ==

The covariance matrix C for the above data values (ignoring the read noise)
would be:

===  ===  ===  ===  ===  ===  ===  ===
--------------------------------------
  1    1    1    1    1    1    1    1
  1    2    2    2    2    2    2    2
  1    2    3    3    3    3    3    3
  1    2    3  104  104  104  104  104
  1    2    3  104  105  105  105  105
  1    2    3  104  105  206  206  206
  1    2    3  104  105  206  207  207
  1    2    3  104  105  206  207  208
===  ===  ===  ===  ===  ===  ===  ===

The read noise was left out of C for simplicity.  To include the read
noise, the square of the read noise divided by the number of frames per
group would be added to each term on the main diagonal.

The solution is the vector P::

    P = (X.T * C^-1 * X)^-1 * (X.T * C^-1 * Y)

where X.T is the transpose of X.  The variances of the parameters P
are given by the diagonal of this matrix::

    (X.T * C^-1 * X)^-1

The GLS solution is computed using three iterations (currently).  For the
first iteration, the actual data are used for populating the covariance
matrix.  For the second and third iterations, the previous fit is used
to populate the covariance matrix.  The covariance matrix should contain
values that actually represent the variance of the data and the covariance
between groups.  If the input data are bad, the covariance matrix might
not contain reasonable estimates of the variances of the data, and this can
result in a very poor fit, along with unrealistic (e.g. negative) variances
for the fitted parameters.  As example of this is when the data are falling
due to severe saturation, but they were not flagged as such because the
data exceeded saturation level before the first read.  The covariance
matrix is constructed under the assumption of an accumulating ramp, and
that is not consistent with rapidly falling data.

The optional output file has a different format when generalized least
squares was specified.  The extension names are YINT, SIGYINT, PEDESTAL,
CRMAG, and SIGCRMAG.  There is one set of values for each integration,
and even if there is only one integration in the exposure, the integration
number is included in the dimensions of the data.  For example, YINT,
SIGYINT, and PEDESTAL have shape (n_int, ny, nx) (which is (nx,ny,n_int)
in IRAF notation), where n_int is the number of integrations in the
exposure, ny is the number of image lines, and nx is the number of image
columns.  YINT is the Y-intercept, the fitted ramp (different for each
pixel) extrapolated back to zero time.  SIGYINT is the error estimate for
YINT.  The time of the first group is taken to be frame_time * (M + 1) / 2,
where frame_time is the time to read out one frame, and M is the number of
frames that were averaged to make a group.  PEDESTAL is the extrapolation
of the first group back to zero time, using the fitted slope.  Note that
PEDESTAL and YINT are similar but not the same.  CRMAG is the magnitude
of each jump (cosmic-ray hit), and SIGCRMAG is the error estimate for
CRMAG.  CRMAG and SIGCRMAG have shape (n_int, ny, nx, max_cr), where
max_cr is the maximum number of jumps (but at least one) that were
identified (by the jump step) for any ramp in any integration in the
exposure.  These arrays are zero-padded, and many elements will likely be
zero.  For a given integration n and pixel [y, x], there will be max_cr
elements for the amplitudes of the jumps.  The amplitude for the first
jump detected in that ramp will be in element [n, y, x, 0], the amplitude
for the second jump will be in element [n, y, x, 1], etc., regardless of
where within the ramp these jumps were detected.  If j jumps were found in
the ramp for pixel [y, x] in integration n, elements [n, y, x, j:max_cr]
will be zero.  The default name for this file will be based on the name
of the input file, using suffix fitoptgls, but the user can override this
name by using the parameter opt_name.

For both ordinary least squares and generalized least squares, pixels
are processed simultaneously in blocks using the array-based functionality
of numpy.  The size of the block depends on the image size and the number
of groups.

Upon successful completion of this step, the status keyword S_RAMP will be
set to COMPLETE.

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
The ramp fitting step has five optional arguments that can be set by the user:

* ``--save_opt``: A True/False value that specifies whether to write
  optional output information.

* ``--opt_name``: A string that can be used to override the default name
  for the optional output information.

* ``--int_name``: A string that can be used to override the default name
  for the integration-by-integration slopes, for the case that the input
  file contains more than one integration.

* ``--algorithm``: A string that can be set to ``GLS`` to mean that
  generalized least squares should be used for the fit.  The default value
  is ``OLS``, which means to use ordinary least squares.

* ``--weighting``: A string that can be set to ``OPTIMAL`` to perform the 
  fitting with the weighting scheme used by Fixsen et al, PASP,112, 1350. 
  This option is only available when using ordinary least squares. The default 
  is ``UNWTD``, which means a uniform weighting scheme will be used.  
