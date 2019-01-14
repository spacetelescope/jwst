.. _calwebb_spec3:

calwebb_spec3: Stage 3 Spectroscopic Processing
===============================================

:Config: calwebb_spec3.cfg
:Class: `~jwst.pipeline.Spec3Pipeline`

Stage 3 processing for spectroscopic observations is intended for combining the 
calibrated data from multiple exposures (e.g. a dither/nod pattern) into a single
combined 2D or 3D spectral product and a combined 1D spectrum.
Before being combined, the exposures may receive additional corrections for the
purpose of background matching and outlier rejection.
The steps applied by the ``calwebb_spec3`` pipeline are shown below.
This pipeline is intended for non-TSO spectra only. TSO spectral data should be
processed using the :ref:`calwebb_tso3 <calwebb_tso3>` pipeline.

.. |c| unicode:: U+2713 .. checkmark

+---------------------------------------------------+-----+-----+-----+-----+-----+--------+--------+
| Instrument/Mode                                   |      NIRSpec    |    MIRI   | NIRISS | NIRCam |
+---------------------------------------------------+-----+-----+-----+-----+-----+--------+--------+
| Step                                              | FS  | MOS | IFU | FS  | MRS | WFSS   | WFSS   |
+===================================================+=====+=====+=====+=====+=====+========+========+
| :ref:`mrs_imatch <mrs_imatch_step>`               |     |     |     |     | |c| |        |        |
+---------------------------------------------------+-----+-----+-----+-----+-----+--------+--------+
| :ref:`outlier_detection <outlier_detection_step>` | |c| | |c| | |c| | |c| | |c| |  |c|   |  |c|   |
+---------------------------------------------------+-----+-----+-----+-----+-----+--------+--------+
| :ref:`resample_spec <resample_step>`              | |c| | |c| |     | |c| |     |  |c|   |  |c|   |
+---------------------------------------------------+-----+-----+-----+-----+-----+--------+--------+
| :ref:`cube_build <cube_build_step>`               |     |     | |c| |     | |c| |        |        |
+---------------------------------------------------+-----+-----+-----+-----+-----+--------+--------+
| :ref:`extract_1d <extract_1d_step>`               | |c| | |c| | |c| | |c| | |c| |  |c|   |  |c|   |
+---------------------------------------------------+-----+-----+-----+-----+-----+--------+--------+

.. note:: The resampling and combining of WFSS spectra is not yet implemented. It
          is expected to be included in DMS Build 7.3 in mid-2019.

Arguments
---------

The ``calwebb_spec3`` pipeline does not have any optional arguments.

Inputs
------

2D calibrated data
^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`,
             `~jwst.datamodels.SlitModel`, or `~jwst.datamodels.MultiSlitModel`
:File suffix: _cal

The inputs to ``calwebb_spec3`` should be in the form of an ASN file that
lists the multiple exposures to be processed into combined output products.
The individual exposures should be calibrated the ("_cal") products from
``calwebb_spec2`` processing.

Outputs
-------

CR-flagged exposures
^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`
:File suffix: _crf

If the :ref:`outlier_detection <outlier_detection_step>` step is applied, a new version of
each input calibrated exposure is created, in which the DQ array has been updated to
flag pixels detected as outliers. These files use the "_crf" (CR-Flagged)
product type suffix and also includes the association candidate ID as a
new field in the original product root name, e.g.
"jw96090001001_03101_00001_nrs2_o001_crf.fits."


2D resampled and combined spectral data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.DrizProductModel`
:File suffix: _s2d

When processing non-IFU modes, a resampled/rectified 2D product of type
"_s2d" is created containing the rectified and combined data for a given
slit/source, which is the output of the :ref:`resample_spec <resample_step>` step.

3D resampled and combined spectral data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.IFUCubeModel`
:File suffix: _s3d

When processing IFU exposures, a resampled/rectified 3D IFU cube product
created by the :ref:`cube_build <cube_build_step>` step is saved as an "_s3d" file.

1D extracted spectral data
^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.MultiSpecModel`
:File suffix: _x1d

All types of inputs result in a 1D extracted spectral data product, which is
saved as a "_x1d" file, and is the result of performing the
:ref:`extract_1d <extract_1d_step>` step on the combined "_s2d" or "_s3d" product.
