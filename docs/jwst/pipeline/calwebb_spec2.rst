.. _calwebb_spec2:
.. _calwebb_tso-spec2:

calwebb_spec2: Stage 2 Spectroscopic Processing
===============================================

:Config: calwebb_spec2.cfg, calwebb_tso-spec2.cfg
:Class: `~jwst.pipeline.Spec2Pipeline`

The ``Spec2Pipeline`` applies additional instrumental corrections and calibrations
to countrate products that result in a fully calibrated individual exposure.
There are two unique configuration files to be used to control this pipeline,
depending on whether the data are to be treated as Time Series Observation (TSO).
Non-TSO exposures use the ``calwebb_spec2`` configuration, which applies all
applicable steps to the data. The ``calwebb_tso-spec2`` configuration, on the other
hand, should be used for TSO exposures, for which some steps are set to be skipped
by default (see the list of steps in the table below). Both configurations call the
``Spec2Pipeline`` module; the only difference is which steps are applied.

The ``Spec2Pipeline`` is the "Swiss army knife" of pipeline modules, containing many
steps that are only applied to certain instruments or instrument modes. The logic for
determining which steps are appropriate is built into the pipeline module itself and
is mostly based on either the instrument name or the exposure type (EXP_TYPE keyword)
of the data.
The list of steps is shown in the table
below and indicates which steps are applied to various spectroscopic modes, as
well as which ones are used by the ``calwebb_tso-spec2.cfg`` configuration for TSO data.
The instrument mode abbreviations used in the table are as follows:

- NIRSpec FS = Fixed Slit
- NIRSpec MOS = Multi-Object Spectroscopy
- NIRSpec IFU = Integral Field Unit
- MIRI FS = LRS Fixed Slit
- MIRI SL = LRS Slitless
- MIRI MRS = Medium Resolution Spectroscopy (IFU)
- NIRISS SOSS = Single Object Slitless Spectroscopy
- NIRISS and NIRCam WFSS = Wide-Field Slitless Spectroscopy

.. |c| unicode:: U+2713 .. checkmark

+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| Instrument/Mode                               |      NIRSpec    |      MIRI       |    NIRISS   | NIRCam | All |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| Step                                          | FS  | MOS | IFU | FS  | SL  | MRS | SOSS | WFSS | WFSS   | TSO |
+===============================================+=====+=====+=====+=====+=====+=====+======+======+========+=====+
| :ref:`assign_wcs <assign_wcs_step>`           | |c| | |c| | |c| | |c| | |c| | |c| |  |c| | |c|  |  |c|   | |c| |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`background <background_step>`           | |c| | |c| | |c| | |c| |     | |c| |  |c| | |c|  |  |c|   |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`imprint <imprint_step>`                 |     | |c| | |c| |     |     |     |      |      |        |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`msaflagopen <msaflagopen_step>`         |     | |c| | |c| |     |     |     |      |      |        |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`extract_2d <extract_2d_step>`\ :sup:`1` | |c| | |c| |     |     |     |     |      | |c|  |  |c|   | |c| |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`srctype <srctype_step>`\ :sup:`1`       | |c| | |c| | |c| | |c| | |c| | |c| |  |c| |      |        | |c| |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`wavecorr <wavecorr_step>`               | |c| | |c| |     |     |     |     |      |      |        |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`flat_field <flatfield_step>`\ :sup:`1`  | |c| | |c| | |c| | |c| | |c| | |c| |  |c| | |c|  |  |c|   | |c| |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`straylight <straylight_step>`           |     |     |     |     |     | |c| |      |      |        |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`fringe <fringe_step>`                   |     |     |     |     |     | |c| |      |      |        |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`pathloss <pathloss_step>`               | |c| | |c| | |c| |     |     |     |  |c| |      |        |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`barshadow <barshadow_step>`             |     | |c| |     |     |     |     |      |      |        |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`photom <photom_step>`                   | |c| | |c| | |c| | |c| | |c| | |c| |  |c| | |c|  |  |c|   | |c| |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`resample_spec <resample_step>`          | |c| | |c| |     | |c| |     |     |      |      |        |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`cube_build <cube_build_step>`           |     |     | |c| |     |     | |c| |      |      |        |     |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`extract_1d <extract_1d_step>`           | |c| | |c| | |c| | |c| | |c| | |c| |  |c| | |c|  |  |c|   | |c| |
+-----------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+

:sup:`1`\ The exact order of the :ref:`extract_2d <extract_2d_step>`, :ref:`srctype <srctype_step>`,
and :ref:`flat_field <flatfield_step>` steps depends on the observing mode.
For NIRISS and NIRCam WFSS, as well as NIRCam TSO grism exposures, the order is
flat_field followed by extract_2d (no wavecorr or srctype).
For all other modes the order is extract_2d, srctype, wavecorr, and flat_field.

The :ref:`resample_spec <resample_step>` step produces a resampled/rectified product for
non-IFU modes of some spectroscopic exposures. If the :ref:`resample_spec <resample_step>` step
is not applied to a given exposure, the :ref:`extract_1d <extract_1d_step>` operation will be
performed on the original (unresampled) data. The :ref:`cube_build <cube_build_step>` step produces
a resampled/rectified cube for IFU exposures, which is then used as input to
the :ref:`extract_1d <extract_1d_step>` step.

Arguments
---------
The ``calwebb_spec2`` pipeline has one optional argument::

  --save_bsub  boolean  default=False

If set to ``True``, the results of the background subtraction step will be saved
to an intermediate file, using a product type of "_bsub" or "_bsubints", depending on
whether the data are 2D (averaged over integrations) or 3D (per-integration results).

Inputs
------

2D or 3D countrate data
^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`,
             or `~jwst.datamodels.CubeModel`
:File suffix: _rate or _rateints

The input to the ``Spec2Pipeline`` pipeline is a countrate exposure, in the form
of either "_rate" or "_rateints" data. A single input file can be processed or an
ASN file listing multiple inputs can be used, in which case the processing steps
will be applied to each input exposure, one at a time.
If "_rateints" products are used as input, each step applies its algorithm to each
integration in the exposure, where appropriate.

Note that the steps :ref:`background <background_step>` and :ref:`imprint <imprint_step>`
can only be executed when the pipeline is given an ASN file as input, because they rely on
multiple, associated exposures to perform their tasks. The ASN file must list not only the
input science exposure(s), but must also list the exposures to be used as background
or imprint.

Background subtraction for Wide-Field Slitless Spectroscopy (WFSS) exposures,
on the other hand, is accomplished by scaling and subtracting a master background
image contained in a CRDS reference file and hence does not require an ASN as input.

The input data model type `~jwst.datamodels.IFUImageModel` is only used for MIRI MRS
and NIRSpec IFU exposures.

Outputs
-------

2D or 3D background-subtracted data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`,
              or `~jwst.datamodels.CubeModel`
:File suffix: _bsub or _bsubints

This is an intermediate product that is only created if "--save_bsub" is set
to ``True`` and will contain the data as output from the :ref:`background <background_step>`
step. If the input is a "_rate" product, this will be a "_bsub" product, while
"_rateints" inputs will be saved as "_bsubints."

2D or 3D calibrated data
^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`,
             `~jwst.datamodels.CubeModel`,
             `~jwst.datamodels.SlitModel`, or `~jwst.datamodels.MultiSlitModel`
:File suffix: _cal or _calints

The output is a fully calibrated, but unrectified, exposure, using the product
type suffix "_cal" or "_calints", dependening on the type of input,
e.g. "jw80600012001_02101_00003_mirimage_cal.fits." This is the output of the
:ref:`photom <photom_step>` step, or whichever step is performed last before applying
either :ref:`resample_spec <resample_step>`, :ref:`cube_build <cube_build_step>`, or
:ref:`extract_1d <extract_1d_step>`.

The output data model type can be any of the 4 listed above and is completely
dependent on the type of input data and the observing mode. For data sets that
do **not** go through :ref:`extract_2d <extract_2d_step>` processing, the output will be
either a `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`, or
`~jwst.datamodels.CubeModel`, matching the corresponding input data type.

Of the data types that do go through :ref:`extract_2d <extract_2d_step>` processing,
the output type will consist of either a single slit model or a mutli-slit model:

- NIRSpec Bright-Object and NIRCam TSO Grism: `~jwst.datamodels.SlitModel`
- NIRSpec Fixed Slit and MOS, as well as WFSS: `~jwst.datamodels.MultiSlitModel`

The multi-slit model is simply an array of multiple slit models, each one
containing the data and relevant meta data for a particular extracted slit or
source. A `~jwst.datamodels.MultiSlitModel` product will contain multiple
tuples of SCI, ERR, DQ, WAVELENGTH, etc. arrays; one for each of the
extracted slits/sources.

2D resampled data
^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.DrizProductModel` or `~jwst.datamodels.MultiProductModel`
:File suffix: _s2d

If the input is a 2D exposure type that gets resampled/rectified by the
:ref:`resample_spec <resample_step>` step, the rectified 2D spectral product is saved as a
"_s2d" file. This image is intended for use as a quick-look product only and is
not used in subsequent processing. The 2D unresampled, calibrated ("_cal")
products are passed along as input to subsequent Stage 3 processing.

If the input to the :ref:`resample_spec <resample_step>` step is a `~jwst.datamodels.MultiSlitModel`,
then the resampled output will be in the form of a
`~jwst.datamodels.MultiProductModel`, which contains an array of individual models,
one per slit. Otherwise the output will be a `~jwst.datamodels.DrizProductModel`.

3D resampled (IFU cube) data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.IFUCubeModel`
:File suffix: _s3d

If the data are NIRSpec IFU or MIRI MRS, the result of the :ref:`cube_build <cube_build_step>`
step will be 3D IFU spectroscopic cube saved to a "_s3d" file. The IFU cube is built from
the data contained in a single exposure and is intended for use as a quick-look
product only. The 2D unresampled, calibrated ("_cal") products are passed along as
input to subsequent Stage 3 processing.

1D extracted spectral data
^^^^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.MultiSpecModel`
:File suffix: _x1d or _x1dints

All types of inputs result in a 1D extracted spectral data product, which is saved
as a "_x1d" or "_x1dints" file, depending on the input type. Observing modes
such as MIRI LRS fixed slit and MRS, NIRCam and NIRISS WFSS, and NIRSpec
fixed slit, MOS, and IFU result in an "_x1d" product containing extracted spectral
data for one or more slits/sources. TSO modes, such as MIRI LRS slitless, NIRCam
TSO grism, NIRISS SOSS, and NIRSpec Bright Object, for which the data are 3D
stacks of integrations, result in "_x1dints" products containing extracted
spectral data for each integration with the exposure.
