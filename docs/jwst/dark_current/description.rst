Description
===========

:Class: `jwst.dark_current.DarkCurrentStep`
:Alias: dark_current

Assumptions
-----------
It is assumed that the input science data have *NOT* had the zero group (or
bias) subtracted. We also do not want the dark subtraction process to remove
the bias signal from the science exposure, therefore the dark reference data
should have their own group zero subtracted from all groups. This means that
group zero of the dark reference data will effectively be zero-valued.

Algorithm
---------
The algorithm for this step is called from the external package ``stcal``, an STScI
effort to unify common calibration processing algorithms for use by multiple observatories.

The dark current step removes dark current from an exposure by subtracting
dark current data stored in a dark reference file in CRDS.

The current implementation uses dark reference files that have been
constructed from exposures using NFRAMES=1 and GROUPGAP=0 (i.e. one
frame per group and no dropped frames) and the maximum number of frames
allowed for an integration. If the science exposure that's being processed
also used NFRAMES=1 and GROUPGAP=0, then the dark reference file data
are directly subtracted group-by-group from the science exposure.

If the science exposure used NFRAMES>1 or GROUPGAP>0, the dark
reference file data are reconstructed on-the-fly by the step to match the frame
averaging and groupgap settings of the science exposure. The reconstructed dark
data are created by averaging NFRAMES adjacent dark frames and skipping
GROUPGAP intervening frames.

The frame-averaged dark is constructed using the following scheme:

* SCI arrays are computed as the mean of the original dark SCI arrays
* ERR arrays are computed as the uncertainty in the mean, using
  :math:`\frac{\sqrt {\sum \mathrm{ERR}^2}}{nframes}`

The dark reference data are not integration-dependent for most instruments,
hence the same group-by-group dark current data are subtracted from every
integration of the science exposure. An exception to this rule is the JWST
MIRI instrument, for which the dark signal **is** integration-dependent, at
least to a certain extent. MIRI dark reference file data is therefore
4-dimensional (ncols x nrows x ngroups x nintegrations). Typical MIRI dark
reference files contain data for only 2 or 3 integrations, which are directly
subtracted from the corresponding first few integrations of the science exposure.
The data in the last integration of the dark reference file is applied to all
remaining science integrations.

The ERR arrays of the science data are currently not modified by this step.

The DQ flags from the dark reference file are propagated into the science
exposure PIXELDQ array using a bitwise OR operation.

Upon successful completion of the dark subtraction the S_DARK keyword is
set to "COMPLETE".

Special Handling
++++++++++++++++
Any pixel values in the dark reference data that are set to NaN will have their
values reset to zero before being subtracted from the science data, which
will effectively skip the dark subtraction operation for those pixels.

**Note**: If the input science exposure contains more groups than the available
dark reference file, no dark subtraction will be applied and the input data
will be returned unchanged.

Subarrays
---------
It is assumed that dark current will be subarray-dependent, therefore this
step makes no attempt to extract subarrays from the dark reference file to
match input subarrays. It instead relies on the presence of matching subarray
dark reference files in CRDS.

JWST/NIRCam Target Acq Subarrays
--------------------------------
Due to the very large number of available NIRCam target acquisition (TA) subarrays,
the instrument team has chosen to not provide dark reference files for any of
the TA subarrays in CRDS.
Requests from the calibration pipeline to CRDS for matching dark reference
files to use when processing a NIRCam TA will result in a reference file name of
"N/A" being returned, which causes the dark subtraction step to skip processing.
Hence dark current will not be subtracted from NIRCam TA subarray exposures.
