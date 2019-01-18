.. _calwebb_image2:
.. _calwebb_tso-image2:

calwebb_image2: Stage 2 Imaging Processing
==========================================

:Config: calwebb_image2.cfg, calwebb_tso-image2.cfg
:Class: `~jwst.pipeline.Image2Pipeline`

Stage 2 imaging processing applies additional instrumental corrections and
calibrations that result in a fully calibrated individual exposure. There are two
unique configuration files to be used to control this pipeline, depending on whether
the data are to be treated as Time Series Observation (TSO). Non-TSO exposures use
the ``calwebb_image2`` configuration, which applies all applicable steps to the data.
The ``calwebb_tso-image2.cfg`` configuration, on the other hand, should be used for
TSO exposures, for which some steps are set to be skipped by default (see the list of
steps in the table below). Both configurations call the ``Image2Pipeline``; the only
difference is which steps are applied.

The list of steps applied by the ``Image2Pipeline`` pipeline is shown in the table
below and indicates which steps are used by the ``calwebb_tso-image2.cfg`` configuration.

.. |check| unicode:: U+2713 .. checkmark

+--------------------------------------+--------------------+
| calwebb_image2                       | calwebb_tso-image2 |
+======================================+====================+
| :ref:`background <background_step>`  |                    |
+--------------------------------------+--------------------+
| :ref:`assign_wcs <assign_wcs_step>`  | |check|            |
+--------------------------------------+--------------------+
| :ref:`flat_field <flatfield_step>`   | |check|            |
+--------------------------------------+--------------------+
| :ref:`photom <photom_step>`          | |check|            |
+--------------------------------------+--------------------+
| :ref:`resample <resample_step>` [1]_ |                    |
+--------------------------------------+--------------------+

.. [1] Resampling is only performed for exposure types "MIR_IMAGE", "NRC_IMAGE", and
   "NIS_IMAGE".

Arguments
---------

The ``calwebb_image2`` pipeline has one optional argument::

  --save_bsub  boolean  default=False

If set to ``True``, the results of
the background subtraction step will be saved to an intermediate file,
using a product type of "_bsub" or "_bsubints", depending on whether the
data are 2D (averaged over integrations) or 3D (per-integration results).

Inputs
------

2D or 3D countrate data
^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
:File suffix: _rate or _rateints

The input to ``Image2Pipeline`` is
a countrate exposure, in the form of either "_rate" or "_rateints"
data. A single input file can be processed or an ASN file listing
multiple inputs can be used, in which case the processing steps will be
applied to each input exposure, one at a time. If "_rateints" products are
used as input, each step applies its algorithm to each
integration in the exposure, where appropriate.

TSO and coronagraphic exposures are expected to use 3D data as input, to be
processed on a per-integration basis. TSO exposures should use the
``calwebb_tso-image2.cfg`` configuration, while coronagraphic exposures
should use the regular ``calwebb_image2.cfg`` configuration.

Outputs
-------

2D or 3D background-subtracted data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
:File suffix: _bsub or _bsubints

This is an intermediate product that is only created if "--save_bsub" is set
to ``True`` and will contain the data as output from the
:ref:`background <background_step>` step.
If the input is a "_rate" product, this will be a "_bsub" product, while
"_rateints" inputs will be saved as "_bsubints."

2D or 3D calibrated data
^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
:File suffix: _cal or _calints

The output is a fully calibrated, but unrectified, exposure, using
the product type suffix "_cal" or "_calints", depending on the type of
input, e.g. "jw80600012001_02101_00003_mirimage_cal.fits".

2D resampled image
^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.DrizProductModel`
:File suffix: _i2d

This is the output of the :ref:`resample <resample_step>` step and is only created
for regular direct imaging observations (not for TSO or coronagraphy 3D data sets).
The output file will use the "_i2d" product type suffix, e.g.
"jw80600012001_02101_00003_mirimage_i2d.fits". Note that this product is
intended for quick-look use only and is not passed along as input to Stage 3
processing. Calibrated, but unrectified (_cal) products are used as input to
Stage 3.

