.. _calwebb_wfs-image2:

calwebb_wfs-image2: Stage 2 WFS&C Processing
============================================

:Config: calwebb_wfs-image2.cfg
:Class: `~jwst.pipeline.Image2Pipeline`

Stage 2 processing of Wavefront Sensing and Control (WFS&C) images duplicates the
processing applied to regular science imaging, with the exception of image resampling.
The ``calwebb_wfs-image2.cfg`` configuration utilizes the regular ``Image2Pipeline``
module, with the :ref:`resample <resample_step>` step set to be skipped, because the
analyis of WFS&C data must be done in the original unrectified image space.
The list of steps is shown in the table below.

.. |check| unicode:: U+2713 .. checkmark

+--------------------------------------+--------------------+
| calwebb_image2                       | calwebb_wfs-image2 |
+======================================+====================+
| :ref:`background <background_step>`  | |check|            |
+--------------------------------------+--------------------+
| :ref:`assign_wcs <assign_wcs_step>`  | |check|            |
+--------------------------------------+--------------------+
| :ref:`flat_field <flatfield_step>`   | |check|            |
+--------------------------------------+--------------------+
| :ref:`photom <photom_step>`          | |check|            |
+--------------------------------------+--------------------+
| :ref:`resample <resample_step>`      |                    |
+--------------------------------------+--------------------+

Arguments
---------
The ``calwebb_wfs-image2`` pipeline does not have any optional arguments.

Inputs
------

2D countrate data
^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _rate

The input to ``Image2Pipeline`` is a countrate exposure, in the form of "_rate" 
data. A single input file can be processed or an ASN file listing
multiple inputs can be used, in which case the processing steps will be
applied to each input exposure, one at a time.

Outputs
-------

2D calibrated data
^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _cal

The output is a fully calibrated, but unrectified, exposure, using
the product type suffix "_cal."
