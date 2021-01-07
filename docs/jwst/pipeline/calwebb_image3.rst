.. _calwebb_image3:

calwebb_image3: Stage 3 Imaging Processing
==========================================

:Config: calwebb_image3.cfg
:Class: `~jwst.pipeline.Image3Pipeline`

Stage 3 processing for direct-imaging observations is intended for combining the 
calibrated data from multiple exposures (e.g. a dither or mosaic pattern) into a
single rectified (distortion corrected) product.
Before being combined, the exposures receive additional corrections for the
purpose of astrometric alignment, background matching, and outlier rejection.
The steps applied by the ``calwebb_image3`` pipeline are shown below.
This pipeline is intended for non-TSO imaging only. TSO imaging data should be
processed using the :ref:`calwebb_tso3 <calwebb_tso3>` pipeline.

+--------------------------------------------------+
| calwebb_image3                                   |
+==================================================+
| :ref:`tweakreg <tweakreg_step>`                  |
+--------------------------------------------------+
| :ref:`skymatch <skymatch_step>`                  |
+--------------------------------------------------+
| :ref:`outlier_detection <outlier_detection_step>`|
+--------------------------------------------------+
| :ref:`resample <resample_step>`                  |
+--------------------------------------------------+
| :ref:`source_catalog <source_catalog_step>`      |
+--------------------------------------------------+

Arguments
---------

The ``calwebb_image3`` pipeline does not have any optional arguments.

Inputs
------

2D calibrated images
^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _cal

The inputs to the ``calwebb_image3`` pipeline are one or more
:ref:`calwebb_image2 <calwebb_image2>`
calibrated ("_cal") image products. In order to process and combine multiple
images, an ASN file must be used as input, listing the exposures to be
processed. It is also possible use a single "_cal" file as input to
``calwebb_image3``, in which case only the :ref:`resample <resample_step>` and
:ref:`source_catalog <source_catalog_step>` steps will be applied.

Outputs
-------

CR-flagged exposures
^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _crf

If the :ref:`outlier_detection <outlier_detection_step>` step is applied, a new version
of each input calibrated exposure is created, in which the DQ array has been updated to
flag pixels detected as outliers. These files use the "_crf" (CR-Flagged)
product type suffix and also includes the association candidate ID as a
new field in the original product root name, e.g.
"jw96090001001_03101_00001_nrca2_o001_crf.fits."

Resampled and combined 2D image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.DrizProductModel`
:File suffix: _i2d

A resampled 2D image product of type "_i2d" is created containing the
combined, rectified association of exposures, which is the direct output of
the :ref:`resample <resample_step>` step.

Source catalog
^^^^^^^^^^^^^^

:Data model: N/A
:File suffix: _cat

The source catalog produced by the :ref:`source_catalog <source_catalog_step>` step
from the "_i2d" product is saved as an ASCII file in ``ecsv`` format, with a product type
of "_cat."
