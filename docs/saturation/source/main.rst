Description
============

The ``saturation`` step flags saturated pixel values. It loops over all
integrations within an exposure, examining them group-by-group, comparing the
science exposure values with defined saturation thresholds for each pixel.
When it finds a pixel value in a given group that is above the threshold, it
sets the ``SATURATED`` flag in the corresponding location of the GROUPDQ
array in the science exposure.

Reference Files
===============
This step requires a SATURATION reference file, which is used to specify the
saturation threshold for each pixel. The saturation files are FITS format,
with 2 IMAGE extensions: ``SCI`` and ``DQ``. They are both 2-D integer arrays.
The values in the SCI array give the saturation threshold in units of DN for
each pixel. The saturation reference file also contains a ``DQ_DEF`` table
extension, which lists the bit assignments for the flag conditions used in
the DQ array.

For pixels having a saturation threshold set to NaN in the reference file,
those thresholds will be replaced by 100000, a very high value that exceeds
any possible science data pixel value. This ensures that these pixels will
not be flagged by this step as saturated. The associated groupdq values will 
be flagged as NO_SAT_CHECK in the step output. Similarly, for pixels flagged 
as NO_SAT_CHECK in the reference file, they will be added to the dq mask, and 
have their saturation values set to be so high they will not be flagged as 
saturated.

The saturation reference files are selected based on instrument, detector and,
where necessary, subarray.

Subarrays
=========
The step will accept either full-frame or subarray saturation reference files.
If only a full-frame reference file is available, the step will extract
subarrays to match those of the science exposure. Otherwise, subarray-specific
saturation reference files will be used if they are available.
