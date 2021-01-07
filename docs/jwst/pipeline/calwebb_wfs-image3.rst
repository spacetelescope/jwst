.. _calwebb_wfs-image3:

calwebb_wfs-image3: Stage 3 WFS&C Processing
============================================

:Config: calwebb_wfs-image3.cfg
:Class: `~jwst.wfs_combine.WfsCombineStep`

Stage 3 processing of Wavefront Sensing and Control (WFS&C) images is only performed
for dithered pairs of WFS&C exposures. The processing applied is not truly a
"pipeline", but consists only of the single :ref:`wfs_combine <wfs_combine_step>` step.
The ``calwebb_wfs-image3.cfg`` configuration exists only for consistency and
compatibility with stage 3 processing of other observing modes. The same result could
be obtained by just running the :ref:`wfs_combine <wfs_combine_step>` step directly.

Arguments
---------
The ``calwebb_wfs-image3`` pipeline has one optional argument::

  --do_refine  boolean  default=False

If set to ``True``, offsets between the dithered images computed from the WCS will be
refined emperically using a cross-correlation technique.
See :ref:`wfs_combine <wfs_combine_step>` for details.

Inputs
------

2D calibrated images
^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _cal

The input to ``calwebb_wfs-image3`` is a pair of calibrated ("_cal") exposures, specified
via an ASN file.

Outputs
-------

2D combined image
^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _wfscmb

The output is a combined image, using the product type suffix "_wfscmb."
See :ref:`wfs_combine <wfs_combine_step>` for details on how this combined
image is produced.
