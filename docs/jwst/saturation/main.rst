Description
============

The ``saturation`` step flags saturated pixel values. It loops over all
integrations within an exposure, examining each one group-by-group, comparing the
science exposure values with defined saturation thresholds for each pixel.
When it finds a pixel value in a given group that is above the threshold, it
sets the ``SATURATED`` flag in the corresponding location of the ``GROUPDQ``
array in the science exposure. It also flags all subsequent groups for that
pixel as saturated. For example, if there are 10 groups in an integration and
group 7 is the first one to cross the saturation threshold for a given pixel,
then groups 7 through 10 will all be flagged for that pixel.

For pixels having a saturation threshold set to NaN in the reference file,
the threshold will be set to 100000, a very high value that exceeds
any possible science data pixel value. This ensures that these pixels will
not be flagged by this step. These pixels will be flagged with
NO_SAT_CHECK in the output ``GROUPDQ`` array. Similarly, pixels flagged
with NO_SAT_CHECK in the reference file ``DQ`` array will have saturation
checking skipped and will have the NO_SAT_CHECK flag set in the output
``PIXELDQ`` array.

Subarrays
=========
The step will accept either full-frame or subarray saturation reference files.
If only a full-frame reference file is available, the step will extract a
subarray to match that of the science exposure. Otherwise, subarray-specific
saturation reference files will be used if they are available.
