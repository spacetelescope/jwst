Description
===========

:Class: `jwst.badpix_selfcal.BadpixSelfcalStep`
:Alias: badpix_selfcal

Overview
--------
The ``badpix_selfcal`` step flags bad pixels in the input data using a self-calibration 
technique based on median filtering along the spectral axis. 
When background exposures are available, those are used to identify bad pixels; 
when unavailable, the science data itself is used.
This correction is applied to `~jwst.datamodels.IFUImageModel` data
directly after the :ref:`assign_wcs <assign_wcs_step>` correction has been applied
in the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline.

Input details
-------------
The input data must be in the form of a `~jwst.datamodels.IFUImageModel` or 
a `~jwst.datamodels.ModelContainer` containing exactly one
science exposure and any number of background exposures. A fits or association file 
that can be read into one of these data models is also acceptable.

Algorithm
---------
The algorithm relies on the assumption that bad pixels are outliers in the data along
the spectral axis. The algorithm proceeds as follows:
* A combined background image is created. If multiple background exposures are available, 
the pixelwise minimum of all background exposures is taken. If only one background exposure
is available, it is used without modification. If no background exposures are 
available, the science data itself is passed in without modification, serving as the 
"background image" for the rest of the procedure, i.e., true self-calibration.
* The combined background image is median-filtered, ignoring NaNs, along the spectral (Y-) axis 
with a user-specified kernel size. The default kernel size is 15 pixels.
* The difference between the original background image and the median-filtered background image
is taken. The highest- and lowest-flux pixels in this difference image are
flagged as bad pixels. The default fraction of pixels to flag is 0.1% of the total number of pixels
on each of the high-flux and low-flux ends of the distribution. This fraction can be adjusted
using the ``flagfrac`` parameter. The total fraction of flagged pixels is thus 2x ``flagfrac``.
* The bad pixels are flagged in the input data by setting the DQ flag to "WARM".
* The bad pixels are also flagged in each background exposure, if available.

Output product
--------------
The output is a new copy of the input `~jwst.datamodels.IFUImageModel`, with the
bad pixels flagged. If the optional ``save_flagged_bkgd`` parameter is set to
``True``, the background exposures with the flagged pixels will also be saved
as separate files. If the entire ``calwebb_spec2`` pipeline is run, the background
exposures passed into the ``background`` step will include the flags from the
``badpix_selfcal`` step.
