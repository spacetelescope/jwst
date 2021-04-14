.. _calwebb_guider:

calwebb_guider: Guide Star Processing
=====================================
:Class: `jwst.pipeline.GuiderPipeline`
:Alias: calwebb_guider

The guider pipeline is only for use with data resulting from the FGS guiding functions:
Identification (ID), Acquisition (ACQ1 and ACQ2), Track, and Fine Guide.
The pipeline applies three detector-level correction and calibration steps to uncalibrated
guider data, as listed in the table below.

+-------------------------------------+
| calwebb_guider                      |
+=====================================+
| :ref:`dq_init <dq_init_step>`       |
+-------------------------------------+
| :ref:`guider_cds <guider_cds_step>` |
+-------------------------------------+
| :ref:`flat_field <flatfield_step>`  |
+-------------------------------------+

Arguments
---------
The ``calwebb_guider`` pipeline does not have any optional arguments.

Inputs
------

4D raw data
^^^^^^^^^^^

:Data model: `~jwst.datamodels.GuiderRawModel`
:File suffix: _uncal

The input to ``calwebb_guider`` is a single raw guide-mode data file, which contains
the original raw data from all of the detector readouts performed during the guider
mode episode. The organization of the 4D data array is analogous to that of 4D raw
science data, having dimensions of ncols x nrows x ngroups x nintegrations.

Outputs
-------

3D calibrated data
^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.GuiderCalModel`
:File suffix: _cal

The output is a 3D (ncols x nrows x nints) countrate product that has been flat-fielded
and has bad pixels flagged. For ID mode data, there is only 1
countrate image produced by the :ref:`guider_cds <guider_cds_step>` step, therefore the
length of the 3rd array axis is 1. For all other modes, there will be a stack of 
multiple countrate images, one per integration. See the :ref:`guider_cds <guider_cds_step>`
step information for details on how the countrate images are produced for each mode.

