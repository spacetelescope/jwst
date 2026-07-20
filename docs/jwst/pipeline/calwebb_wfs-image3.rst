.. _calwebb_wfs-image3:

calwebb_wfs-image3: Stage 3 WFS&C Processing
============================================

:Class: `jwst.wfs_combine.wfs_combine_step.WfsCombineStep`
:Alias: calwebb_wfs-image3

Stage 3 processing of Wavefront Sensing and Control (WFS&C) images is only performed
for dithered pairs of WFS&C exposures. The processing applied is not truly a
"pipeline", but consists only of the single :ref:`wfs_combine <wfs_combine_step>` step.
The ``calwebb_wfs-image3`` alias exists only for consistency and
compatibility with stage 3 processing of other observing modes. The same result could
be obtained by just running the :ref:`wfs_combine <wfs_combine_step>` step directly.

For more details:

* :ref:`wfs_combine_args`
* :ref:`wfs_comb_inputs`
* :ref:`wfs_comb_outputs`
* :ref:`wfs_comb_algo`
