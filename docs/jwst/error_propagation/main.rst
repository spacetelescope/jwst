Description
===========
 
Steps in the pipelines calculate arrays for variances, or modify variance arrays 
from previous steps in the pipeline.  In some cases these arrays are only used 
internally to the current step.  For several steps, these arrays must be propagated to
subsequent steps in the pipeline, or included in the cumulate error array.

Stage 1 pipelines 
-----------------
Stage 1 pipelines perform detector-level corrections and ramp fitting for
individual exposures, for nearly all imaging and spectroscopic modes. Details 
of the pipelines can be found at `Stage 1 Pipelines
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_detector1.html>`_

There are two unique configuration files to be used to control this pipeline, 
depending on whether the data are to be treated as a Time Series Observation 
(in which case some of the steps are set to be skipped by default).

The steps to be executed differ for Near-IR and MIRI exposures. For 
Near-IR, the steps are group_scale, dq_init, saturation, superbias, refpix, 
linearity, persistence, dark_current, jump, ramp fitting, and gain scale.
For MIRI, the steps are group_scale, dq_init, saturation, firstframe, lastframe,
linearity, rcsd, dark_current, refpix, jump, ramp fitting, and gain scale.

Of all of those steps, below are the particular steps that alter the ERR,
VAR_POISSON, or VAR_RNOISE arrays of the science ramp data and output those
altered arrays in the science product that is passed to the next step in the
pipeline.  For the steps not listed below, those arrays do not exist or are
propagated unchanged to the next pipeline step.


Ramp_fitting
++++++++++++
This step calculates and populates the VAR_POISSON, VAR_RNOISE arrays
(the variances and errors associated with ramp fitting) in the 'rate' and
'rateints' files, and updates the ERR arrays in those files. The details of the 
calculations can be found at `Ramp fitting 
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/ramp_fitting/main.html>`_



Group_scale
+++++++++++
The group_scale step is applied after the ramp_fitting step, and applies to both the 
rate and rateints products from ramp_fit. The correction is applied to the ERR, 
(and will also be applied to the) VAR_POISSON, and VAR_RNOISE arrays.  The details
of the calculations can be found at `Group Scale
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/group_scale/description.html>`_


Stage 2 pipelines 
-----------------
Stage 2 pipelines perform additional instrument-level and observing-mode corrections and 
calibrations to produce fully calibrated exposures. There are two Stage 2 pipelines: one 
for `imaging 
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_image2.html>`_
and one for `spectroscopy 
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_spec2.html>`_
In these pipelines, the various steps that apply corrections and calibrations
apply those same corrections/calibrations to all variance arrays arrays.

There are two unique configuration files to be used to control these pipelines, 
depending on whether the data are to be treated as a Time Series Observation 
(in which case some of the steps are set to be skipped by default).


Stage 2 Imaging Processing pipeline - calwebb_image2 & calwebb_tso-image2
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The steps executed in this pipeline depend on whether the data are to be treated
as TSO.  For non-TSO exposures, the steps are background image subtraction,
assign_wcs, flat_field, photom, and resample.  For TSO exposures, the steps are
assign_wcs, flat_field, and photom.

Of the steps listed above, below are the steps that alter the ERR or VAR arrays
of the science ramp data and output those altered arrays in the science product 
that is passed to the next step in the pipeline.  For the steps not listed below,
those arrays do not exist or are propagated unchanged to the next pipeline step.

Flat_field
~~~~~~~~~~
(NOTE: get from the latest online readthedocs for this step once James updates it.)
The SCI array from the flat-field reference file is divided into the VAR_POISSON
and VAR_RNOISE arrays of the science data set. The variance in the flat is
calculated in this step; the ERR array is then recalculated as these three 
variances added in quadrature.  The newly created VAR_FLAT array and the updated 
VAR_POISSON and VAR_RNOISE are output to the science products(s). For details of
the flat fielding correction, see `Flatfielding
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/index.html>`_


Photom
~~~~~~
The photom step is applied for the modes: imaging and non-IFU spectroscopy,
NIRSpec IFU, and MIRI MRS. The application of the step to these modes is very
similar - for details see `Photom
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/photom/index.html>`_

For these modes, the calibration information includes a scalar
flux conversion constant, as well as optional arrays of wavelength and
relative response (as a function of wavelength).  The combination of the scalar
conversion factor and the 2-D response values is applied to the science data,
including the SCI and ERR arrays, as well as the variance (VAR_POISSON,
VAR_RNOISE, and VAR_FLAT) arrays. The correction values are divided into the SCI 
and ERR arrays, and the square of the correction values are divided into the 
variance arrays.


Stage 2 Spectroscopic Processing pipeline - calwebb_spec2 & calwebb_tso-spec2
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The steps executed in this pipeline depend on whether the data are to be treated
as TSO.  For non-TSO exposures, the steps are assign_wcs, background image
subtraction, imprint, msaflagopen, extract_2d, flat_field, srctype, straylight,
fringe, pathloss, barshadow, photom, resample_spec, cube_build, and extract_1d.
For TSO exposures, some of the steps are in the list above are skipped.  The
paricular steps executed depend on particular instruments or instrument modes-
see `Calwebb_spec2
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_spec2.html#calwebb-spec2>`_

Of the steps listed above, below are the steps that alter the ERR or VAR
arrays of the science ramp data and output those altered arrays in the science
product that is passed to the next step in the pipeline.  For the steps not
listed below, those arrays do not exist or are propagated unchanged to the next
pipeline step.

Flat_field
~~~~~~~~~~
{NOTE get from the latest online readthedocs for the step once James updates it.}
The SCI array from the flat-field reference file is divided into the VAR_POISSON
and VAR_RNOISE arrays of the science data set. The variance in the flat is
calculated in this step; the ERR array is then recalculated as these three
variances added in quadrature.  The newly created VAR_FLAT array and the updated
VAR_POISSON and VAR_RNOISE are output to the science products(s). For details of
the flat fielding correction, see `Flatfielding
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/index.html>`_


Fringe 
~~~~~~
For MIRI MRS (IFU) mode exposures, the ERR array in the science data set is
divided by a fringe reference image.  For details of the fringe correction, see 
`Fringe
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/fringe/index.html>`_


Barshadow 
~~~~~~~~~
This step is applied only to NIRSpec MSA data for extended sources. Once the
2-D correction array for each slit has been computed, it is applied to the
science (SCI), error (ERR), and variance (VAR_POISSON, VAR_RNOISE, and VAR_FLAT)
data arrays of the slit.  The correction values are divided into the SCI and ERR
arrays, and the square of the correction values are divided into the variance 
arrays.   For details of the bar shadow correction, see `Barshadow
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/barshadow/index.html>`_


Photom
~~~~~~ 
The calibration information includes a scalar flux conversion constant, as well as
optional arrays of wavelength and relative response (as a function of wavelength).
The combination of the scalar conversion factor and the 2-D response values is
applied to the science data, including the SCI and ERR arrays, as well as the
variance (VAR_POISSON, VAR_RNOISE, and VAR_FLAT) arrays. The correction values
are divided into the SCI and ERR arrays, and the square of the correction values 
are divided into the variance arrays.  For details of the photom correction, see
`Photom
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/photom/index.html>`_


Cube build
~~~~~~~~~~
In the output spectral cube, the ERR extension contains the uncertainty on the 
SCI values.  For details of the cube build step, see `CubeBuild
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/cube_build/index.html>`_


Stage 2 WFS&C Processing pipeline - calwebb_wfs-image2
++++++++++++++++++++++++++++++++++++++++++++++++++++++
This pipeline processes Wavefront Sensing and Control (WFS&C) images, which
duplicates the processing applied to regular science imaging, with the exception
of image resampling. The steps in this pipeline are background image subtraction, 
assign_wcs, flat_field, photom, and resample.  Of these steps, below are the steps 
that alter the ERR or VAR arrays of the science ramp data and output those altered
arrays in the science product that is passed to the next step in the pipeline.
For the steps not listed below, those arrays do not exist or are propagated
unchanged to the next pipeline step.


Flat_field
~~~~~~~~~~
{NOTE get from the latest online readthedocs for the step once James updates it.}
The SCI array from the flat-field reference file is divided into the VAR_POISSON
and VAR_RNOISE arrays of the science data set. The variance in the flat is
calculated in this step; the ERR array is then recalculated as these three
variances added in quadrature.  The newly created VAR_FLAT array and the updated
VAR_POISSON and VAR_RNOISE are output to the science products(s).  For details of
the flat fielding correction, see `Flatfielding
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/index.html>`_


Photom
~~~~~~ 
The calibration information includes a scalar flux conversion constant, as well as
optional arrays of wavelength and relative response (as a function of wavelength).
The combination of the scalar conversion factor and the 2-D response values is
applied to the science data, including the SCI and ERR arrays, as well as the
variance (VAR_POISSON, VAR_RNOISE, and VAR_FLAT) arrays. The correction values
are divided into the SCI and ERR arrays, and the square of the correction values 
are divided into the variance arrays.  For details of the photom correction, see 
`Photom
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/photom/index.html>`_


Stage 3 pipelines
-----------------
Stage 3 pipelines perform corrections that work with multiple exposures and in
most cases produce some kind of combined product).  There are unique pipelines
for imaging, spectroscopic, coronographic, AMI and TSO observations. OMIT ??>> For details
of these pipelines see <add link here>.


Stage 3 Image Processing Pipeline - calwebb_image3
++++++++++++++++++++++++++++++++++++++++++++++++++
This pipeline is for non-TSO imaging only, and combines the calibrated data 
from multiple exposures (e.g. a dither or mosaic pattern) into a single rectified
(distortion corrected) product.  The steps in this pipeline are tweakreg, 
skymatch, outlier_detection, resample, source_catalog.  The only one of these 
steps that either uses or modifies variance/error arrays that are propagated
through the pipeline is outlier_detection, described here:

Outlier_detection 
~~~~~~~~~~~~~~~~~
This step uses the ERR array to make a local noise model, based on the readnoise
and calibration errors of earlier steps in the pipeline. This step does not modify
the ERR array or any of the VAR arrays.


Stage 3 Spectroscopic Processing Pipeline - calwebb_spec3
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This pipeline is for non-TSO imaging only, and is intended for combining the
calibrated data from multiple exposures (e.g. a dither/nod pattern) into a
single combined 2D or 3D spectral product and a combined 1D spectrum.  The steps 
in this pipeline are master_background, exp_to_source, mrs_match, 
outlier_detection, resample_spec, cube_build, extract_1d, and combine_1d.  The
paricular steps executed depend on particular instruments or instrument modes; 
for details see `Calwebb_spec3
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_spec3.html>`_

Of the steps listed above, below are the steps that either use or modify
variance/error arrays that are propagated through these pipelines.

Outlier_detection 
~~~~~~~~~~~~~~~~~
This step uses the ERR array to make a local noise model, based on the readnoise
and calibration errors of earlier steps in the pipeline. This step does not modify
the ERR array or any of the VAR arrays.

Cube_build
~~~~~~~~~~
Cube_build takes MIRI or NIRSpec IFU calibrated 2-D images and produces 3-D
spectral cubes.  In the output spectral cube, the SCI exension contains the
surface brightness of cube spaxels in units of mJy/arcsecond^2, and the ERR
extension contains the uncertainty on the SCI values.


Stage 3 Coronographic Processing Pipeline - calwebb_coron3 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The stage 3 coronagraphic pipeline is to be applied to associations of
calibrated NIRCam coronagraphic and MIRI Lyot and 4QPM exposures, and is used 
to produce PSF-subtracted, resampled, combined images of the source object.  The
steps in this pipeline are stack_refs, align_refs, klip, outlier_detection, and
resample.  The only one of these steps that either uses or modifies variance/error
arrays that are propagated through the pipeline is outlier_detection, described below.
For more details on this pipeline see `Calwebb_coron3
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_coron3.html>`_


Outlier_detection
~~~~~~~~~~~~~~~~~
This step uses the ERR array to make a local noise model, based on the readnoise
and calibration errors of earlier steps in the pipeline. This step does not modify
the ERR array or any of the VAR arrays. 


Stage 3 Time-Series Observation (TSO) Processing Pipeline - calwebb_tso3
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The stage 3 TSO pipeline is to be applied to associations of calibrated TSO
exposures.  The steps in this pipeline are outlier_detection, tso_photometry,
extract_1d, and white_light.  The only one of these steps that either uses or
modifies variance/error arrays that are propagated through the pipeline is
outlier_detection, described here:

Outlier_detection 
~~~~~~~~~~~~~~~~~
This step uses the ERR array to make a local noise model, based on the readnoise
and calibration errors of earlier steps in the pipeline. This step does not modify
the ERR array or any of the VAR arrays.

Stage 3 WFS&C Processing pipeline - calwebb_wfs-image3
++++++++++++++++++++++++++++++++++++++++++++++++++++++
Stage 3 processing of Wavefront Sensing and Control (WFS&C) images is only
performed for dithered pairs of WFS&C exposures. The processing applied is not 
truly a “pipeline”, but consists only of the single wfs_combine step. This step
could modify variance/error arrays, but only if the optional 'do_refine' is set
to True (which is *NOT* the default in pipeline use). In this case the "refined 
algorithm" will be used, and the ERR array values will be altered on output.


Other pipelines
---------------
In addtion to the Level 1, 2, and 3 pipelines discussed so far, there are 2 other 
pipelines- Dark Processing, and Guide Star Processing.

Dark Processing pipeline - calwebb_dark
+++++++++++++++++++++++++++++++++++++++
The Dark Pipeline applies basic detector-level corrections to all dark exposures.
For details on thie pipeline see `Calwebb_dark
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_dark.html>`_
All of the steps in this pipeline are group_scale, dq_init, saturation, ipc, superbias, 
refpix, linearity, and rcsd.  The only one of these steps that either uses or
modifies variance/error arrays that are propagated through the pipeline is
group_scale, described here:

Group_scale
~~~~~~~~~~~
This step is applied after the ramp_fitting step, and applies to both the rate and
rateints products from ramp_fit. The correction is applied to the ERR, and will
also be applied to the VAR_POISSON, and VAR_RNOISE arrays.  The details of the 
calculations can be found at `Group Scale
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/group_scale/description.html>`_


Guide Star Processing pipeline - calwebb_guider
+++++++++++++++++++++++++++++++++++++++++++++++
The Guider pipeline is only for use with data resulting from the FGS guiding functions.
The steps in this pipeline are dq_init, guider_cds, and flat_field.
 
The only one of these steps that either uses or modifies variance/error arrays
that are propagated through the pipeline is flat_field, described here:

Flat_field
~~~~~~~~~~
{NOTE get from the latest online readthedocs for the step once James updates it.}
The SCI array from the flat-field reference file is divided into the VAR_POISSON
and VAR_RNOISE arrays of the science data set. The variance in the flat is
calculated in this step; the ERR array is then recalculated as these three
variances added in quadrature.  The newly created VAR_FLAT array and the updated
VAR_POISSON and VAR_RNOISE are output to the science products(s).  For details of
the flat fielding correction, see `Flatfielding
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/index.html>`_


The table below is a summary of which steps create and output variance arrays,
modify and output the cumulate error or variance arrays, use locally but do not output arrays,
and which level pipelines each step is in.


================= ======================= ====================================== ============================ ======================
STEP              Creates arrays          Modifies arrays                        Step-specific use of arrays  Pipeline Level
================= ======================= ====================================== ============================ ======================
photom            None                    ERR, VAR_POISSON, VAR_RNOISE, VAR_FLAT None                         Stages 2,3
flat field        VAR_FLAT                ERR, VAR_POISSON, VAR_RNOISE           None                         Stage 2, Guide Star
wfs_combine       None                    ERR                                    None                         Stage 3
group_scale       None                    ERR, VAR_POISSON, VAR_RNOISE           None                         Stage 1, Dark Pipeline
fringe            None                    ERR                                    None                         Stage 2
cube_build        None                    ERR                                    None                         Stages 2,3
bar shadow        None                    ERR, VAR_POISSON, VAR_RNOISE, VAR_FLAT None                         Stage 2
outlier detection None                    None                                   ERR                          Stage 3
ramp_fitting      VAR_POISSON, VAR_RNOISE None                                   None                         Stage 1

================= ======================= ====================================== ============================ ======================



