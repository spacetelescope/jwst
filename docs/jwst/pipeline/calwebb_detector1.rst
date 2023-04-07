.. _calwebb_detector1:
.. _calwebb_tso1:

calwebb_detector1: Stage 1 Detector Processing
==============================================

:Class: `jwst.pipeline.Detector1Pipeline`
:Alias: calwebb_detector1

The ``Detector1Pipeline`` applies basic detector-level corrections to all exposure
types (imaging, spectroscopic, coronagraphic, etc.). It is applied to one
exposure at a time.
It is sometimes referred to as "ramps-to-slopes" processing, because the input raw data
are in the form of one or more ramps (integrations) containing accumulating
counts from the non-destructive detector readouts and the output is a corrected
countrate (slope) image.

There are two general configurations for this pipeline, depending on whether the
data are to be treated as a Time Series Observation (TSO). The configuration is
provided by CRDS and the reftype ``pars-detector1pipeline``. In general, for
Non-TSO exposures, all applicable steps are applied to the data. For TSO
exposures, some steps are set to be skipped by default (see the list of steps in
the table below).

The list of steps applied by the ``Detector1Pipeline`` pipeline is shown in the
table below. Note that MIRI exposures use some instrument-specific steps and
some of the steps are applied in a different order than for Near-IR (NIR)
instrument exposures.

Several steps in this pipeline include special handling for NIRCam "Frame 0"
data. The NIRCam instrument has the ability to downlink the image from the
initial readout that follows the detector reset at the start of each integration
in an exposure. These images are distinct from the first group of each integration
when on-board frame averaging is done. In these cases, the first group contains
data from multiple frames, while frame zero is always composed of just the
first frame following the reset. It can be used to recover an estimated slope for
pixels that go into saturation already in the first group (see more details on
that process in the :ref:`ramp_fitting <ramp_fitting_step>` step). In order for
the frame zero image to be utilized during ramp fitting, it must have all of
the same calibrations and corrections applied as the first group in the various
``Detector1Pipeline`` steps. This includes the :ref:`saturation <saturation_step>`,
:ref:`superbias <superbias_step>`, :ref:`refpix <refpix_step>`, and
:ref:`linearity <linearity_step>` steps. Other steps do not have a direct effect
on either the first group or frame zero pixel values.

.. |check| unicode:: U+2713 .. checkmark

+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| Near-IR                                                                                  | MIRI                                                        |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| Step                                                                 | Non-TSO | TSO     | Step                                    | Non-TSO | TSO     |
+======================================================================+=========+=========+=========================================+=========+=========+
| :ref:`group_scale <group_scale_step>`                                | |check| | |check| | :ref:`group_scale <group_scale_step>`   | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`dq_init <dq_init_step>`                                        | |check| | |check| | :ref:`dq_init <dq_init_step>`           | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`saturation <saturation_step>`                                  | |check| | |check| | :ref:`saturation <saturation_step>`     | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`ipc <ipc_step>` [1]_                                           |         |         | :ref:`ipc <ipc_step>`                   |         |         |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`superbias <superbias_step>`                                    | |check| | |check| | :ref:`firstframe <firstframe_step>`     | |check| |         |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`refpix <refpix_step>`                                          | |check| | |check| | :ref:`lastframe <lastframe_step>`       | |check| |         |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
|                                                                      |         |         | :ref:`reset <reset_step>`               | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`linearity <linearity_step>`                                    | |check| | |check| | :ref:`linearity <linearity_step>`       | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`persistence <persistence_step>` [2]_                           | |check| |         | :ref:`rscd <rscd_step>`                 | |check| |         |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`dark_current <dark_current_step>`                              | |check| | |check| | :ref:`dark_current <dark_current_step>` | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
|                                                                      |         |         | :ref:`refpix <refpix_step>`             | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+----------------------+------------------+---------+---------+
| :ref:`jump <jump_step>`                                              | |check| | |check| | :ref:`jump <jump_step>`                 | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`undersampling_correction <undersampling_correction_step>` [3]_ | |check| |         |                                         |         |         |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`ramp_fitting <ramp_fitting_step>`                              | |check| | |check| | :ref:`ramp_fitting <ramp_fitting_step>` | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+
| :ref:`gain_scale <gain_scale_step>`                                  | |check| | |check| | :ref:`gain_scale <gain_scale_step>`     | |check| | |check| |
+----------------------------------------------------------------------+---------+---------+-----------------------------------------+---------+---------+


.. [1] By default, the parameter reference `pars-detector1pipeline`
   retrieved from CRDS will skip the :ref:`ipc <ipc_step>` step for all instruments.
.. [2] The :ref:`persistence <persistence_step>` step is currently hardwired to be skipped in
   the `Detector1Pipeline` module for all NIRSpec exposures.
.. [3] By default, the :ref:`undersampling_correction <undersampling_correction_step>` step is skipped in
   the `Detector1Pipeline` module for all instruments.
   

Arguments
---------
The ``calwebb_detector1`` pipeline has one optional argument::

  --save_calibrated_ramp  boolean  default=False

If set to ``True``, the pipeline will save intermediate data to a file as it
exists at the end of the :ref:`jump <jump_step>` step. The data
at this stage of the pipeline are still in the form of the original 4D ramps
(ncols x nrows x ngroups x nints) and have had all of the detector-level
correction steps applied to it, including the detection and flagging of
Cosmic-Ray (CR) hits within each ramp (integration). If created, the name of the
intermediate file will be constructed from the root name of the input file, with
the new product type suffix "_ramp" appended,
e.g. "jw80600012001_02101_00003_mirimage_ramp.fits".

Inputs
------

4D raw data
+++++++++++

:Data model: `~jwst.datamodels.RampModel`
:File suffix: _uncal

The input to ``Detector1Pipeline`` is a single raw exposure,
e.g. "jw80600012001_02101_00003_mirimage_uncal.fits", which contains the
original raw data from all of the detector readouts in the exposure
(ncols x nrows x ngroups x nintegrations).

Note that in the operational environment, the
input will be in the form of a `~jwst.datamodels.Level1bModel`, which only
contains the 4D array of detector pixel values, along with some optional
extensions. When such a file is loaded into the pipeline, it is immediately
converted into a `~jwst.datamodels.RampModel`, and has all additional data arrays
for errors and Data Quality flags created and initialized to zero.

The input can also contain a 3D cube of NIRCam "Frame 0" data, where
each image plane in the 3D cube is the initial frame for each integration in
the exposure. Only present when the option to downlink the frame zero data
was selected in the observing program.

Outputs
-------

4D corrected ramp
+++++++++++++++++

:Data model: `~jwst.datamodels.RampModel`
:File suffix: _ramp

Result of applying all pipeline steps up through the :ref:`jump <jump_step>` step,
to produce corrected and CR-flagged 4D ramp data, which will have the same data dimensions
as the input raw 4D data (ncols x nrows x ngroups x nints). Only created when the
pipeline argument ``--save_calibrated_ramp`` is set to ``True`` (default is ``False``).

2D countrate product
++++++++++++++++++++

:Data model: `~jwst.datamodels.ImageModel` or `~jwst.datamodels.IFUImageModel`
:File suffix: _rate

All types of inputs result in a 2D countrate product,
based on averaging over all of the integrations within the exposure.
The output file will be of type "_rate", e.g.
"jw80600012001_02101_00003_mirimage_rate.fits". The 2D "_rate" product is passed along
to subsequent pipeline modules for all non-TSO and non-Coronagraphic exposures.
For MIRI MRS and NIRSpec IFU exposures, the output data model will be
`~jwst.datamodels.IFUImageModel`, while all others will be `~jwst.datamodels.ImageModel`.

3D countrate product
++++++++++++++++++++

:Data model: `~jwst.datamodels.CubeModel`
:File suffix: _rateints

A 3D countrate product is created that contains the individual
results of each integration. The 2D countrate images for each integration are
stacked along the 3rd axis of the data cubes (ncols x nrows x nints). This
output file will be of type "_rateints". The 3D "_rateints" product is passed along
to subsequent pipeline modules for all TSO and Coronagraphic exposures.

.. include:: ../references_general/pars-detector1pipeline_reffile.inc
