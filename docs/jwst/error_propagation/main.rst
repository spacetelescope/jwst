Description
-----------
 
Steps in the pipelines calculate arrays for variances, or modify variance arrays 
from previous steps in the pipeline.  In some cases these arrays are only used 
internally to the current step.  For several steps, these arrays must be propagated to
subsequent steps in the pipeline, or included in the cumulative error array. The
cumulative error (ERR) array values are always recomputed as the quadratic sum
of all variances whenever a new type of variance array is added to the pipeline
flow or whenever any of the variance arrays is updated by a step.


Stage 1 pipelines 
-----------------
Stage 1 pipelines perform detector-level corrections and ramp fitting for
individual exposures, for nearly all imaging and spectroscopic modes. Details 
of the pipelines can be found at :ref:`Stage 1 Pipelines <calwebb_detector1>`.

Of all of the Stage 1 pipeline steps, below are the particular steps that alter the ERR,
VAR_POISSON, or VAR_RNOISE arrays of the science countrate data and output those
altered arrays in the science product that is passed to the next step in the
pipeline.  For the steps not listed below, those arrays do not exist or are
propagated unchanged to the next pipeline step.


Ramp_fitting
++++++++++++
This step calculates and populates the VAR_POISSON and VAR_RNOISE arrays
(the variances associated with ramp fitting) in the 'rate' and
'rateints' files, and updates the ERR arrays in those files. The details of the 
calculations can be found at :ref:`ramp_fitting <ramp_fitting_step>`.


Gain_scale
++++++++++
The gain_scale step is applied after the ramp_fitting step, and applies to both the 
rate and rateints products from ramp_fit. The correction is applied to the ERR, 
VAR_POISSON, and VAR_RNOISE arrays.  The details of the calculations can be
found at :ref:`gain_scale <gain_scale_step>`.


Stage 2 pipelines 
-----------------
Stage 2 pipelines perform additional instrument-level and observing-mode corrections and 
calibrations to produce fully calibrated exposures. There are two Stage 2 pipelines: one 
for imaging :ref:`calwebb_image2 <calwebb_image2>` and one for 
spectroscopy :ref:`calwebb_spec2 <calwebb_spec2>`.

In these pipelines, the various steps that apply corrections and calibrations
apply those same corrections/calibrations to all variance arrays.


Stage 2 Imaging Processing pipeline - calwebb_image2, calwebb_tso-image2, calwebb_wfs-image2
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For the Stage 2 imaging pipelines, below are the steps that alter the ERR or
VAR arrays of the science countrate data and output those altered arrays in the
science product that is passed to the next step in the pipeline.  For the steps
not listed below, those arrays do not exist or are propagated unchanged to the
next pipeline step.

Flat_field
~~~~~~~~~~
The SCI array from the flat-field reference file is divided into the VAR_POISSON
and VAR_RNOISE arrays of the science data set. The variance in the flat is
calculated in this step; the ERR array is then recalculated as these three 
variances added in quadrature.  The newly created VAR_FLAT array and the updated 
VAR_POISSON and VAR_RNOISE are output to the science products(s). For details of
the flat fielding correction, see :ref:`flat_field <flatfield_step>`.


Photom
~~~~~~
The flux calibration that's applied to the science data includes a scalar
flux conversion constant, as well as optional arrays of wavelength and
relative response (as a function of wavelength).  The combined calibrations
are applied to the science data, including the SCI and ERR arrays, as well as
the variance (VAR_POISSON, VAR_RNOISE, and VAR_FLAT) arrays. The correction
values are divided into the SCI and ERR arrays, and the square of the correction
values are divided into the variance arrays.  For more details see 
:ref:`photom <photom_step>`.


Stage 2 Spectroscopic Processing pipeline - calwebb_spec2 & calwebb_tso-spec2
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For the Stage 2 spectroscopic pipeline, below are the steps that alter the ERR
or VAR arrays of the science countrate data and output those altered arrays in
the science product that is passed to the next step in the pipeline.  For the
steps not listed below, those arrays do not exist or are propagated unchanged
to the next pipeline step.

Flat_field
~~~~~~~~~~
The SCI array from the flat-field reference file is divided into the VAR_POISSON
and VAR_RNOISE arrays of the science data set. The variance in the flat is
calculated in this step; the ERR array is then recalculated as these three
variances added in quadrature.  The newly created VAR_FLAT array and the updated
VAR_POISSON and VAR_RNOISE are output to the science products(s). For details of
the flat fielding correction, see :ref:`flat_field <flatfield_step>`.


Fringe 
~~~~~~
For MIRI MRS (IFU) mode exposures, the ERR array in the science data set is
divided by a fringe reference image.  For details of the fringe correction, see 
:ref:`fringe <fringe_step>`.


Barshadow 
~~~~~~~~~
This step is applied only to NIRSpec MSA data for extended sources. Once the
2-D correction array for each slit has been computed, it is applied to the
science (SCI), error (ERR), and variance (VAR_POISSON, VAR_RNOISE, and VAR_FLAT)
data arrays of the slit.  The correction values are divided into the SCI and ERR
arrays, and the square of the correction values are divided into the variance 
arrays.   For details of the bar shadow correction, see
:ref:`barshadow <barshadow_step>`.


Photom
~~~~~~ 
The calibration information includes a scalar flux conversion constant, as well as
optional arrays of wavelength and relative response (as a function of wavelength).
The combination of the scalar conversion factor and the 2-D response values is
applied to the science data, including the SCI and ERR arrays, as well as the
variance (VAR_POISSON, VAR_RNOISE, and VAR_FLAT) arrays. The correction values
are divided into the SCI and ERR arrays, and the square of the correction values 
are divided into the variance arrays.  For details of the photom correction, see
:ref:`photom <photom_step>`.


Cube build
~~~~~~~~~~
In the output spectral cube, the ERR extension contains the uncertainty on the 
SCI values.  For details of the cube build step, see
:ref:`cube_build <cube_build_step>`.


Stage 3 pipelines
-----------------
Stage 3 pipelines perform operations that work with multiple exposures and in
most cases produce some kind of combined product.  The operations in these
pipelines that either use or modify variance/error arrays that are propagated 
through the pipeline are outlier_detection, cube_build, and wfs_combine.

Outlier_detection
+++++++++++++++++
The outlier detection step is used in all Stage 3 pipelines.  It uses the ERR array to
make a local noise model, based on the readnoise and calibration errors of earlier 
steps in the pipeline. This step does not modify the ERR array or any of the VAR
arrays.

Cube_build
++++++++++
The cube_build step is only applied in the Stage 3 Spectroscopic (calwebb_spec3)
pipeline.  It takes MIRI or NIRSpec IFU calibrated 2-D images and produces 3-D
spectral cubes.  In the output spectral cube, the SCI exension contains the
surface brightness of cube spaxels in units of mJy/arcsecond^2, and the ERR
extension contains the uncertainty on the SCI values.

Wfs_combine
+++++++++++
The wfs_combine step is only applied in the Stage 3 Wavefront Sensing and Control
(calwebb_wfs-image3) pipeline for dithered pairs of WFS&C exposures.  The processing
applied is not truly a “pipeline”, but consists only of the single wfs_combine
step. This step could modify variance/error arrays, but only if the optional
'do_refine' is set to True (which is *NOT* the default in pipeline use). In this
case the "refined algorithm" will be used, and the ERR array values will be
altered on output.


The table below is a summary of which steps create and output variance arrays,
modify and output the cumulative error or variance arrays, use locally but do not
output arrays, and which level pipeline(s) each step is in. 

================= ======================= ====================================== ============================ =================
STEP              Creates arrays          Modifies arrays                        Step-specific use of arrays  Pipeline Level(s)
================= ======================= ====================================== ============================ =================
ramp_fitting      VAR_POISSON, VAR_RNOISE None                                   None                         Stage 1
gain_scale        None                    ERR, VAR_POISSON, VAR_RNOISE           None                         Stage 1
flat field        VAR_FLAT                ERR, VAR_POISSON, VAR_RNOISE           None                         Stage 2
fringe            None                    ERR                                    None                         Stage 2
bar shadow        None                    ERR, VAR_POISSON, VAR_RNOISE, VAR_FLAT None                         Stage 2
photom            None                    ERR, VAR_POISSON, VAR_RNOISE, VAR_FLAT None                         Stages 2,3
cube_build        None                    ERR                                    None                         Stages 2,3
outlier detection None                    None                                   ERR                          Stage 3
wfs_combine       None                    ERR                                    None                         Stage 3

================= ======================= ====================================== ============================ =================



