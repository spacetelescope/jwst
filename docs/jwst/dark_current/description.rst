Description
===========

Assumptions
-----------

It is assumed that the input science data have *NOT* had the zero group (or
bias) subtracted. We also do not want the dark subtraction process to remove
the bias signal from the science exposure, therefore the dark reference data
should have their own group zero subtracted from all groups. This means that
group zero of the dark reference data will effectively be zero-valued.

Algorithm
---------

The dark current step removes dark current from a JWST exposure by subtracting
dark current data stored in a dark reference file.

The current implementation uses dark reference files that have been
constructed from exposures using ``nframes=1`` and ``groupgap=0`` (i.e. one
frame per group and no dropped frames) and the maximum number of frames
allowed for an integration. If the science exposure that's being processed
also used ``nframes=1`` and ``groupgap=0``, then the dark reference file data
are directly subtracted frame-by-frame from the science exposure.

If the science exposure used ``nframes>1`` or ``groupgap>0``, the dark
reference file data are reconstructed internally to match the frame averaging
and groupgap settings of the science exposure. The reconstructed dark data are
constructed by averaging ``nframes`` adjacent dark frames and skipping
``groupgap`` intervening frames.

The frame-averaged dark is constructed using the following scheme:

* SCI arrays are computed as the mean of the original dark SCI arrays
* ERR arrays are computed as the uncertainty of the mean, using
  :math:`\frac{\sqrt {\sum \mathrm{ERR}^2}}{nframes}`

For each integration in the input science exposure, the averaged dark data are
then subtracted, group-by-group, from the science exposure groups, as follows:

* Each SCI group of the dark data are subtracted from the corresponding SCI
  group of the science data
* The ERR arrays of the science data are not modified

Any pixel values in the dark reference data that are set to NaN will have their
values reset to zero before being subtracted from the science data, which
will effectively skip the dark subtraction operation for those pixels.

The dark DQ array is combined with the science exposure PIXELDQ array using a
bitwise OR operation.

**Note**: If the input science exposure contains more frames than the available
dark reference file, no dark subtraction will be applied and the input data
will be returned unchanged.

Subarrays
---------

It is assumed that dark current will be subarray-dependent, therefore this
step makes no attempt to extract subarrays from the dark reference file to
match input subarrays. It instead relies on the presence of matching subarray
dark reference files in CRDS.
