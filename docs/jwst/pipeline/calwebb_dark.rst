.. _calwebb_dark:

calwebb_dark: Dark Processing
=============================

:Config: calwebb_dark.cfg
:Class: `~jwst.pipeline.DarkPipeline`

The ``DarkPipeline`` applies basic detector-level corrections to all dark exposures.
It is identical to the :ref:`calwebb_detector1 <calwebb_detector1>` pipeline, except
that it stops processing immediately before the :ref:`dark_current <dark_current_step>` step.
The list of steps is shown below. As with the :ref:`calwebb_detector1 <calwebb_detector1>`
pipeline, the order of steps is a bit different for MIRI exposures.

+---------------------------------------+-----------------------------------------+
| Near-IR                               | MIRI                                    |
+=======================================+=========================================+
| :ref:`group_scale <group_scale_step>` | :ref:`group_scale <group_scale_step>`   |
+---------------------------------------+-----------------------------------------+
| :ref:`dq_init <dq_init_step>`         | :ref:`dq_init <dq_init_step>`           |
+---------------------------------------+-----------------------------------------+
| :ref:`saturation <saturation_step>`   | :ref:`saturation <saturation_step>`     |
+---------------------------------------+-----------------------------------------+
| :ref:`ipc <ipc_step>` [1]_            | :ref:`ipc <ipc_step>`                   |
+---------------------------------------+-----------------------------------------+
| :ref:`superbias <superbias_step>`     | :ref:`firstframe <firstframe_step>`     |
+---------------------------------------+-----------------------------------------+
| :ref:`refpix <refpix_step>`           | :ref:`lastframe <lastframe_step>`       |
+---------------------------------------+-----------------------------------------+
| :ref:`linearity <linearity_step>`     | :ref:`linearity <linearity_step>`       |
+---------------------------------------+-----------------------------------------+
|                                       | :ref:`rscd <rscd_step>`                 |
+---------------------------------------+-----------------------------------------+

.. [1] The :ref:`ipc <ipc_step>` step is currently set to be skipped by default in the
   "calwebb_dark.cfg" configuration file for all instruments.

Arguments
---------
The ``calwebb_dark`` pipeline has no optional arguments.

Inputs
------

4D raw data
+++++++++++

:Data model: `~jwst.datamodels.RampModel`
:File suffix: _uncal

The input to ``DarkPipeline`` is a single raw dark exposure,
which contains the original raw data from all of the detector readouts in the exposure
(ncols x nrows x ngroups x nintegrations).

Outputs
-------

4D corrected ramp
+++++++++++++++++

:Data model: `~jwst.datamodels.RampModel`
:File suffix: _dark

Result of applying all pipeline steps listed above.
Will have the same data dimensions as the
input raw 4D data (ncols x nints x ngroups x nints).
