Description
============

This step determines the mean count rate for each pixel by performing a linear
fit to the data in the input (jump) file.  The fit is done using "ordinary
least squares" (the "generalized least squares" is no longer an option).
The fit is performed independently for each pixel.  There are up to three
output files. The primary output file, giving the slope at each pixel, is
always produced.  If the input exposure contains more than one integration, the
resulting slope images from each integration are stored as a data cube in a
second output data product.  A third, optional output product is also available
and is produced only when the step parameter 'save_opt' is True (the default is
False).  The output values will be in units of counts per second.  Following a
description of the fitting algorithm, these three type of output files are
detailed below.


The count rate for each pixel is determined by a linear fit to the
cosmic-ray-free and saturation-free ramp intervals for each pixel; hereafter
this interval will be referred to as a "segment". The fitting algorithm does an 
'optimal' linear fit, which is the weighting used by Fixsen et 
al, PASP,112, 1350. ('unweighted' in which all groups for a pixel are equally
weighted, is no longer a weighting option.) Segments are derived using
the 4-D GROUPDQ array of the input data set, under the assumption that the jump
step will have already flagged CR's. Segments are also terminated where
saturation flags are found. Pixels are processed simultaneously in blocks 
using the array-based functionality of numpy.  The size of the block depends
on the image size and the number of groups.


Each of the next 4 paragraphs describes how special cases are handled.


If the input dataset has only a single group in each integration, the count rate
for all unsaturated pixels in that integration will be calculated to be the
value of the science data in that group divided by the group time.  If the
input dataset has only two groups per integration, the count rate for all
unsaturated pixels in each integration will be calculated using the differences 
of the 2 valid values of the science data.


For datasets having more than a single group in each integration, a ramp having 
a segment with only a single group is processed differently depending on the 
number and size of the other segments in the ramp. If a ramp has only a single
segment having only a single group, the count rate will be calculated to be the
value of the science data in that group divided by the group time.  If a ramp
has a segment having a single group, and at least one other segment having more
than one good group, only data from the segment(s) having more than a single 
good group will be used to calculate the count rate.


The data is checked for ramps in which there is good data in the first group, 
but all first differences for the ramp are undefined because there are too few
good groups past the first. For such ramps the first differences will be set to
equal the data in the first group.  The first difference is used to estimate
the slope of the ramp, as explained in the 'segment-specific computations' 
section below.


If any input dataset contains ramps saturated in their second group, the count
rates for those pixels in that integration will be calculated to be the value
of the science data in the first group divided by the group time. 


For all input datasets, including those special cases just described,
arrays for the primary output product are computed as follows.


After computing the slopes for all segments for a given pixel, the final slope is
determined as a weighted average from all segments in all integrations, and is
written to a file as the primary output product.  In this output product, the
4-D GROUPDQ from all integrations is compressed into 2-D, which is then merged
(using a bitwise OR) with the input 2-D PIXELDQ to create the output DQ array. 
The 3-D VAR_POISSON and VAR_RNOISE arrays from all integrations are averaged
into corresponding 2-D output arrays.  If the ramp fitting step is run by
itself, the output file name will have the suffix '_RampFit' or the suffix
'_RampFitStep'; if the ramp fitting step is run as part of the calwebb_detector1
pipeline, the final output file name will have the suffix '_rate'.  In either
case, the user can override this name by specifying an output file name.


If the input exposure contains more than one integration, the resulting slope
images from each integration are stored as a data cube in a second output data
product.  Each plane of the 3-D SCI, ERR, DQ, VAR_POISSON, and VAR_RNOISE
arrays in this product is the result for a given integration.  In this output
product, the GROUPDQ data for a given integration is compressed into 2-D, which
is then merged with the input 2-D PIXELDQ to create the output DQ array for each
integration. The 3-D VAR_POISSON and VAR_RNOISE from an integration are
calculated by averaging over the fit segments in the corresponding 4-D arrays.
By default, the name of this output product is based on the name of the input
file and will have the suffix '_rateints'; the user can override this name by
specifying a name using the parameter int_name.


A third, optional output product is also available and is produced only when
the step parameter 'save_opt' is True (the default is False).  This optional
product contains 4-D arrays called SLOPE, SIGSLOPE, YINT, SIGYINT, WEIGHTS,
VAR_POISSON, and VAR_RNOISE which contain the slopes, uncertainties in the
slopes, y-intercept, uncertainty in the y-intercept, fitting weights, the
variance of the slope due to poisson noise only, and the variance of the slope
due to read noise only for each segment of each pixel. (Calculation of the two
variance arrays requires retrieving readnoise and gain values from their
respective reference files.)  The y-intercept refers to the result of the fit
at an exposure time of zero.  This product also contains a 3-D array called
PEDESTAL, which gives the signal at zero exposure time for each pixel, and the
4-D CRMAG array, which contains the magnitude of each group that was flagged as
having a CR hit.  By default, the name of this output file is based on the name
of the input file and will have the suffix '_fitopt'; the user can override 
this name by specifying a name using the parameter opt_name.  In this optional
output product, the pedestal array is calculated for each integration by
extrapolating the final slope (the weighted average of the slopes of all of
ramp segments in the integration) for each pixel from its value at the first
group to an exposure time of zero. Any pixel that is saturated on the first
group is given a pedestal value of 0. Before compression, the cosmic ray
magnitude array is equivalent to the input SCI array but with the only nonzero
values being those whose pixel locations are flagged in the input GROUPDQ as
cosmic ray hits. The array is compressed, removing all groups in which all the
values are 0 for pixels having at least one group with a non-zero magnitude.
The order of the cosmic rays within the ramp is preserved.


Slopes and their variances are calculated for each segment, for each integration,
and for the entire dataset. As defined above, a segment is a set of contiguous
groups where none of the groups are saturated or cosmic ray-affected.  The 
appropriate slopes and variances are output to the primary output product, the 
integration-specific output product, and the optional output product. The 
following is a description of these computations. The notation in the equations
is the following: the type of noise (when appropriate) will appear as the
superscript 'R', 'P', or 'C' for readnoise, Poisson noise, or combined,
respectively; and the form of the data will appear as the subscript: 's', 'i',
'o' for segment, integration, or overall (for the entire dataset), respectively.


Segment-specific computations:
------------------------------

The slope of each segment is calculated using the least-squares method with 
optimal weighting. The variance of the slope of the segment due to read noise is: 

.. math::  
   var^R_{s} = \frac{12 \ R^2 }{ (ngroups_{s}^3 - ngroups_{s})(tgroup^2) } \,,

where :math:`R` is the noise in the difference between 2 frames, 
:math:`ngroups_{s}` is the number of groups in the segment, and :math:`tgroup` is the group 
time in seconds (from the keyword TGROUP).  

The variance of the slope of the segment due to Poisson noise is: 

.. math::  
   var^P_{s} = \frac{ slope_{est} }{  tgroup \times gain\ (ngroups_{s} -1)}  \,,


where :math:`gain` is the gain for the pixel (from the GAIN reference file),
in e/DN. The :math:`slope_{est}` is an overall estimated slope of the pixel,
calculated by taking the median of the first differences of the groups that are
unaffected by saturation and cosmic rays, in all integrations. This is a more
robust estimate of the slope than the segment-specific slope, which may be noisy
for short segments. 

The combined variance of the slope of the segment is the sum of the variances: 

.. math::  
   var^C_{s} = var^R_{s} + var^P_{s}


Integration-specific computations:
----------------------------------  
The variance of the slope for the integration due to read noise is:

.. math::  
   var^R_{i} = \frac{1}{ \sum_{s} \frac{1}{ var^R_{s} }}  \,,

where the sum is over all segments in the integration.

The variance of the slope for the integration due to Poisson noise is: 

.. math::  
   var^P_{i} = \frac{1}{ \sum_{s} \frac{1}{ var^P_{s}}}  

The combined variance of the slope for the integration is due to both Poisson and read
noise: 

.. math::  
   var^C_{i} = \frac{1}{ \sum_{s} \frac{1}{ var^R_{s} + var^P_{s}}}

The slope for the integration depends on the slope and the combined variance of each segment's slope:

.. math::  
   slope_{i} = \frac{ \sum_{s}{ \frac{slope_{s}} {var^C_{s}}}} { \sum_{s}{ \frac{1} {var^C_{s}}}}



Exposure-level computations:
----------------------------

The variance of the slope due to read noise depends on a sum over all integrations: 

.. math::  
   var^R_{o} = \frac{1}{ \sum_{i} \frac{1}{ var^R_{i}}} 

The variance of the slope due to Poisson noise is: 

.. math::  
   var^P_{o} = \frac{1}{ \sum_{i} \frac{1}{ var^P_{i}}}


The combined variance of the slope is the sum of the variances: 

.. math::  
   var^C_{o} = var^R_{o} + var^P_{o}

The square root of the combined variance is what gets stored in the ERR array of
the primary output.



The overall slope depends on the slope and the combined variance of the slope of each integration's segments, so is a sum over integrations and segments:

.. math::    
    slope_{o} = \frac{ \sum_{i,s}{ \frac{slope_{i,s}} {var^C_{i,s}}}} { \sum_{i,s}{ \frac{1} {var^C_{i,s}}}}


Upon successful completion of this step, the status keyword S_RAMP will be set
to COMPLETE.

The MIRI first frame correction step flags all pixels in the first group of data
in each integration of a MIRI exposure having more than 3 groups, so that those 
data do not get used in either the jump detection or ramp fitting steps. 
Similarly, the MIRI last frame correction step flags all pixels in the last 
group of data in each integration of a MIRI exposure having more than 2 groups, 
so that those data do not get used in either the jump detection or ramp fitting 
steps. The ramp fitting will only fit data if there are at least 2 good groups 
of data, and will log a warning otherwise.



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


Error Propagation
=================

Error propagation in the ramp fitting step is implemented by storing the square
root of the exposure-level combined variance in the ERR array of the primary
output product. This combined variance of the exposure-level slope is the sum
of the variance of the slope due to the Poisson noise and the variance of the 
slope due to the read noise. These two variances are also separately written
to the extensions VAR_POISSON and VAR_RNOISE in the primary output.


At the integration-level, the variance of the integration-level slope due to
Poisson noise is written to the VAR_POISSON extension in the
integration-specific product, and the variance of the integration-level slope
due to read noise noise is written to the VAR_RNOISE extension.  The combined
variance of the slope is due to both the Poisson noise and the read noise, and
its square root is written to the ERR extension. 


For the optional output product, the variance of the slope due to the Poisson
noise of the segment-specific slope is written to the VAR_POISSON extension.
Similarly, the variance of the slope due to the read noise of the
segment-specific slope  is written to the VAR_RNOISE extension.
