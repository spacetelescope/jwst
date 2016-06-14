Description
===========

The last frame correction step applies the last frame reference file.
This is a MIRI specific correction. The MIRI arrays are reset sequentially
by row pairs. The last frame of an integration ramp on a given pixel is
influenced by signal coupled through the reset of the adjacent row pair. 
The result is that the odd and even rows both show anomalous offsets in the
last read on an integration. 

The "last frame effect" manifests itself in two ways. One is an effect 
on the bias level, which is seen in all images (darks or with illumination) 
and the other is an effect on any collected signal in the array at the 
time the last frame is being read out. This step only corrects the offset 
in the bias level. 

This step subtracts the lastframe reference data from the last group of a
science data integration. This correction is not integration 
dependent. The DQ flags from the last frame  reference file
are combined with the science PIXELDQ array using numpy's bitwise_or function.
The ERR arrays of the science data are currently not modified at all.

Subarrays
---------

The last frame correction may  be subarray-dependent. The current 
implimentation assumes that it will be so this step makes no attempt to 
extract subarrays from the lastframe reference file to match input subarrays. 
It instead relies on the presence of matching subarray last frame  reference 
files in CRDS. 

