Description
-----------
Steps in the various pipeline modules calculate variances due to different sources of
noise or modify variances that were computed by previous steps.  In some cases the
variance arrays are only used internally within a given step.  For several steps,
these arrays must be propagated to subsequent steps in the pipeline. Anytime a step
creates or updates variances, the total error (ERR) array values are always recomputed
as the square root of the quadratic sum of all variances available at the time.
Note that the ERR array values are always expressed as standard deviation
(i.e. square root of the variance).

The table below is a summary of which steps create or update variance and error arrays,
as well as which steps make use of these data. Details of how each step computes or
uses these data are given in the subsequent sections below.

================= ===== ======================= ====================================== ==========
Step              Stage Creates arrays          Updates arrays                         Step uses
================= ===== ======================= ====================================== ==========
ramp_fitting        1   VAR_POISSON, VAR_RNOISE ERR                                    None
gain_scale          1   None                    ERR, VAR_POISSON, VAR_RNOISE           None
flat_field          2   VAR_FLAT                ERR, VAR_POISSON, VAR_RNOISE           None
fringe              2   None                    ERR                                    None
barshadow           2   None                    ERR, VAR_POISSON, VAR_RNOISE, VAR_FLAT None
pathloss            2   None                    ERR, VAR_POISSON, VAR_RNOISE, VAR_FLAT None
photom              2   None                    ERR, VAR_POISSON, VAR_RNOISE, VAR_FLAT None
outlier_detection   3   None                    None                                   ERR
resample            3   None                    None                                   VAR_RNOISE
wfs_combine         3   None                    ERR                                    None
================= ===== ======================= ====================================== ==========

Stage 1 Pipelines 
-----------------
Stage 1 pipelines perform detector-level corrections and ramp fitting for
individual exposures, for nearly all imaging and spectroscopic modes. Details 
of the pipelines can be found at :ref:`Stage 1 Pipelines <calwebb_detector1>`.

The Stage 1 pipeline steps that alter the ERR, VAR_POISSON, or VAR_RNOISE arrays of
the science countrate data are discussed below.
Any step not listed here does not alter or use the variance or error arrays
in any way and simply propagates the information to the next step.

ramp_fitting
++++++++++++
This step calculates and populates the VAR_POISSON and VAR_RNOISE arrays
in the 'rate' and 'rateints' files, and updates the ERR array as the square root of the
quadratic sum of the variances. VAR_POISSON and VAR_RNOISE represent the uncertainty in the
computed slopes (per pixel) due to Poisson and read noise, respectively.
The details of the calculations can be found at :ref:`ramp_fitting <ramp_fitting_step>`.

gain_scale
++++++++++
The ``gain_scale`` step is applied after ``ramp_fitting``, and applies to both the 
rate and rateints products. The gain correction is applied to the ERR, 
VAR_POISSON, and VAR_RNOISE arrays.  The SCI and ERR arrays are multiplied by the
gain correction factor, and the variance arrays are multiplied by the square of
the gain correction factor. More details can be
found at :ref:`gain_scale <gain_scale_step>`.

Stage 2 Pipelines 
-----------------
Stage 2 pipelines perform additional instrument-level and observing-mode corrections and 
calibrations to produce fully calibrated exposures. There are two main Stage 2 pipelines:
one for imaging :ref:`calwebb_image2 <calwebb_image2>` and one for 
spectroscopy :ref:`calwebb_spec2 <calwebb_spec2>`.
In these pipelines, the various steps that apply corrections and calibrations
apply those same corrections/calibrations to all variance arrays and update the total
ERR array.

flat_field
++++++++++
The SCI array of the exposure being processed is divided by the flat-field reference
image.  The VAR_FLAT array is created, computed from the science data and the flat-field
reference file ERR array. 

For all modes except GUIDER, the VAR_POISSON and VAR_RNOISE arrays are divided by the
square of the flat. The science data ERR array is then updated to be the square root
of the sum of the three variance arrays. 

For the GUIDER mode, their are no VAR_POISSON and VAR_RNOISE arrays. The science data
ERR array is updated to be the square root of the sum of the variance VAR_FLAT and the
square of the incoming science ERR array. 

For more details see :ref:`flat_field <flatfield_step>`.

fringe 
++++++
For MIRI MRS (IFU) mode exposures, the SCI and ERR arrays in the science exposure
are divided by the fringe reference image.  For details of the fringe correction, see 
:ref:`fringe <fringe_step>`.

barshadow 
+++++++++
This step is applied only to NIRSpec MSA data for extended sources. Once the
2-D correction array for each slit has been computed, it is applied to the
science (SCI), error (ERR), and variance (VAR_POISSON, VAR_RNOISE, and VAR_FLAT)
arrays of the slit.  The correction values are divided into the SCI and ERR
arrays, and the square of the correction values are divided into the variance 
arrays.   For details of the bar shadow correction, see
:ref:`barshadow <barshadow_step>`.

pathloss
++++++++
The ``pathloss`` step corrects NIRSpec and NIRISS SOSS data for various types of
light losses. The correction factors are divided into the SCI and ERR arrays of
the science data, and the square of the correction values are divided into the
variance arrays. For details of this step, see :ref:`pathloss <pathloss_step>`.

photom
++++++ 
The calibration information for the ``photom`` step includes a scalar flux conversion
constant, as well as optional arrays of wavelength and relative response (as a
function of wavelength). The combination of the scalar conversion factor and any 2-D
response values is applied to the science data, including the SCI and ERR arrays,
as well as the variance (VAR_POISSON, VAR_RNOISE, and VAR_FLAT) arrays. The flux
calibration values are multiplied into the science exposure SCI and ERR arrays,
and the square of the calibration values is multiplied into all variance arrays.
For details of the photom correction, see :ref:`photom <photom_step>`.

Stage 3 pipelines
-----------------
Stage 3 pipelines perform operations that work with multiple exposures and in
most cases produce some kind of combined product.  The operations in these
pipelines that either use or modify variance/error arrays that are propagated 
through the pipeline are ``outlier_detection`` and ``wfs_combine``.

outlier_detection
+++++++++++++++++
The ``outlier_detection`` step is used in all Stage 3 pipelines.  It uses the ERR array to
make a local noise model, based on the readnoise and calibration errors of earlier 
steps in the pipeline. This step does not modify the ERR array or any of the VAR
arrays.

resample/resample_spec
++++++++++++++++++++++
The ``resample`` and ``resample_spec`` steps make use of the VAR_RNOISE array to
compute weights that are used when combining data with the ``weight_type=ivm``
option selected. The step also resamples all of the variance and error arrays,
using the same output WCS frame as the science data.

wfs_combine
+++++++++++
The ``wfs_combine`` step is only applied in the Stage 3 Wavefront Sensing and Control
(calwebb_wfs-image3) pipeline for dithered pairs of WFS&C exposures.
This step can modify variance/error arrays, but only if the optional
"do_refine" parameter is set to True (the default is False). In this
case the algorithm to refine image offsets will be used and the ERR array values will be
altered on output, using a combination of the input image errors.
See the step documentation at :ref:`wfs_combine <wfs_combine_step>` for
more details.
