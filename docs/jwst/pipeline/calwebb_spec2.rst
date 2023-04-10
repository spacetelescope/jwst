.. _calwebb_spec2:
.. _calwebb_tso-spec2:

calwebb_spec2: Stage 2 Spectroscopic Processing
===============================================

:Class: `jwst.pipeline.Spec2Pipeline`
:Alias: calwebb_spec2

The ``Spec2Pipeline`` applies additional instrumental corrections and
calibrations to countrate products that result in a fully calibrated individual
exposure. There are two general configurations for this pipeline, depending on
whether the data are to be treated as Time Series Observation (TSO). In general,
for non-TSO exposures, all applicable steps are applied to the data. For TSO
exposures, some steps are set to be skipped by default (see the list of steps in
the table below).

The ``Spec2Pipeline`` is the "Swiss army knife" of pipeline modules, containing
many steps that are only applied to certain instruments or instrument modes. The
logic for determining which steps are appropriate is built into the pipeline
module itself and determined by the CRDS ``pars-spec2pipeline`` parameter
reference file. Logic is mostly based on either the instrument name or the
exposure type (EXP_TYPE keyword) of the data.

Science Exposures
-----------------

The list of steps shown in the table below indicates which steps are
applied to various spectroscopic modes for JWST science exposures, including
TSO exposures. The instrument mode abbreviations used in the table are as follows:

- NIRSpec FS = Fixed Slit
- NIRSpec MOS = Multi-Object Spectroscopy
- NIRSpec IFU = Integral Field Unit
- MIRI FS = LRS Fixed Slit
- MIRI SL = LRS Slitless
- MIRI MRS = Medium Resolution Spectroscopy (IFU)
- NIRISS SOSS = Single Object Slitless Spectroscopy
- NIRISS and NIRCam WFSS = Wide-Field Slitless Spectroscopy

.. |c| unicode:: U+2713 .. checkmark

+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| Instrument/Mode                                          |      NIRSpec    |      MIRI       |    NIRISS   | NIRCam | All |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| Step                                                     | FS  | MOS | IFU | FS  | SL  | MRS | SOSS | WFSS | WFSS   | TSO |
+==========================================================+=====+=====+=====+=====+=====+=====+======+======+========+=====+
| :ref:`assign_wcs <assign_wcs_step>`                      | |c| | |c| | |c| | |c| | |c| | |c| |  |c| | |c|  |  |c|   | |c| |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`imprint <imprint_step>`                            |     | |c| | |c| |     |     |     |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`background <background_step>`                      | |c| | |c| | |c| | |c| |     | |c| |  |c| | |c|  |  |c|   |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`msaflagopen <msaflagopen_step>`                    |     | |c| | |c| |     |     |     |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`extract_2d <extract_2d_step>`\ :sup:`1`            | |c| | |c| |     |     |     |     |      | |c|  |  |c|   | |c| |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`srctype <srctype_step>`\ :sup:`1`                  | |c| | |c| | |c| | |c| | |c| | |c| |  |c| | |c|  |  |c|   | |c| |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`master_background <master_background_step>`        |     | |c| |     |     |     |     |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`wavecorr <wavecorr_step>`                          | |c| | |c| |     |     |     |     |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`flat_field <flatfield_step>`\ :sup:`1`             | |c| | |c| | |c| | |c| | |c| | |c| |  |c| | |c|  |  |c|   | |c| |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`straylight <straylight_step>`                      |     |     |     |     |     | |c| |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`fringe <fringe_step>`                              |     |     |     |     |     | |c| |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`pathloss <pathloss_step>`                          | |c| | |c| | |c| | |c| |     |     |  |c| |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`barshadow <barshadow_step>`                        |     | |c| |     |     |     |     |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`wfss_contam <wfss_contam_step>`                    |     |     |     |     |     |     |      | |c|  |  |c|   |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`photom <photom_step>`                              | |c| | |c| | |c| | |c| | |c| | |c| |  |c| | |c|  |  |c|   | |c| |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`residual_fringe <residual_fringe_step>` \ :sup:`2` |     |     |     |     |     | |c| |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`resample_spec <resample_step>`                     | |c| | |c| |     | |c| |     |     |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`cube_build <cube_build_step>`                      |     |     | |c| |     |     | |c| |      |      |        |     |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+
| :ref:`extract_1d <extract_1d_step>`                      | |c| | |c| | |c| | |c| | |c| | |c| |  |c| | |c|  |  |c|   | |c| |
+----------------------------------------------------------+-----+-----+-----+-----+-----+-----+------+------+--------+-----+

:sup:`1`\ The exact order of the :ref:`extract_2d <extract_2d_step>`, :ref:`srctype <srctype_step>`,
and :ref:`flat_field <flatfield_step>` steps depends on the observing mode.
For NIRISS and NIRCam WFSS, as well as NIRCam TSO grism exposures, the order is
flat_field, extract_2d, and srctype (no wavecorr).
For all other modes the order is extract_2d, srctype, wavecorr, and flat_field.

:sup:`2`\ By default the :ref:`residual_fringe <residual_fringe_step>` is skipped in the ``calwebb_spec2`` pipeline. 

Notice that NIRSpec MOS is the only mode to receive master background subtraction
in the ``calwebb_spec2`` pipeline. All other spectral modes have master background
subtraction applied in the :ref:`calwebb_spec3 <calwebb_spec3>` pipeline.

The :ref:`resample_spec <resample_step>` step produces a resampled/rectified product for
non-IFU modes of some spectroscopic exposures. If the :ref:`resample_spec <resample_step>` step
is not applied to a given exposure, the :ref:`extract_1d <extract_1d_step>` operation will be
performed on the original (unresampled) data. The :ref:`cube_build <cube_build_step>` step produces
a resampled/rectified cube for IFU exposures, which is then used as input to
the :ref:`extract_1d <extract_1d_step>` step.

NIRSpec Lamp Exposures
----------------------

The ``Spec2Pipeline`` works slightly differently for NIRSpec lamp exposures.
These are identified by the EXP_TYPE values of NRS_LAMP, NRS_AUTOWAVE or
NRS_AUTOFLAT.  Using the EXP_TYPE keyword in this way means that another keyword
is needed to specify whether the data are Fixed Slit, MOS, IFU or Brightobj.
This is the OPMODE keyword, which maps to the ``jwst.datamodel`` attribute
``.meta.instrument.lamp_mode``.  This keyword can take the following values in
exposures that undergo ``Spec2Pipeline`` processing:

- BRIGHTOBJ = Bright Object mode (uses fixed slits)
- FIXEDSLIT = Fixed slit mode
- IFU = Integral Field Unit mode
- MSASPEC = Multi-Object Spectrograph Mode

OPMODE can also take the values of GRATING-ONLY and NONE, but only in some
engineering-only situations, and can take the value of IMAGE for imaging
data.  None of these values will trigger the execution of the ``Spec2Pipeline``.

NIRSpec calibration lamps are identified by the LAMP keyword,
which maps to the ``jwst.datamodel`` attribute ``.meta.instrument.lamp_state``.
The lamps are either line lamps, used for wavelength calibration, or continuum
lamps, which are used for flatfielding.  Each is paired with a specific grating:

+-----------+---------------------------+-------------------+
| Lamp name | Wavelength range (micron) | Used with grating |
+===========+===========================+===================+
| FLAT1     | 1.0 - 1.8                 |   G140M, G140H    |
+-----------+---------------------------+-------------------+
| FLAT2     | 1.7 - 3.0                 |   G235M, G235H    |
+-----------+---------------------------+-------------------+
| FLAT3     | 2.9 - 5.0                 |   G395M, G395H    |
+-----------+---------------------------+-------------------+
| FLAT4     | 0.7 - 1.4                 |   G140M, G140H    |
+-----------+---------------------------+-------------------+
| FLAT5     | 1.0 - 5.0                 |       PRISM       |
+-----------+---------------------------+-------------------+
| LINE1     | 1.0 - 1.8                 |   G140M, G140H    |
+-----------+---------------------------+-------------------+
| LINE2     | 1.7 - 3.0                 |   G235M, G235H    |
+-----------+---------------------------+-------------------+
| LINE3     | 2.9 - 5.0                 |   G395M, G395H    |
+-----------+---------------------------+-------------------+
| LINE4     | 0.6 - 5.0                 |       PRISM       |
+-----------+---------------------------+-------------------+
| REF       | 1.3 - 1.7                 |   G140M, G140H    |
+-----------+---------------------------+-------------------+

The pairing comes because the calibration unit lightpath doesn't pass through
the filter wheel, so each lamp has its own filter identical to those in the
filter wheel.

The list of ``Spec2Pipeline`` steps to be run for NIRSpec lamp exposures is
shown in the table below and indicates which steps are
applied to various spectroscopic modes. The instrument mode
abbreviations used in the table are as follows:

- NIRSpec FS = Fixed Slit (also Brightobj)
- NIRSpec MOS = Multi-Object Spectroscopy
- NIRSpec IFU = Integral Field Unit

+---------------------------------------+------------+--------------+-----------------+--------------+
|    Pipeline Step                      |         NRS_LAMP          |  NRS_AUTOWAVE   | NRS_AUTOFLAT |
+---------------------------------------+------------+--------------+-----------------+              +
|                                       |   LINE     |     FLAT     |                 |  (MOS only)  |
+=======================================+============+==============+=================+==============+
| :ref:`assign_wcs <assign_wcs_step>`   |   ALL      |     ALL      |       ALL       |       ALL    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`imprint <imprint_step>`         |  NONE      |     IFU      |      NONE       |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`background <background_step>`   |  NONE      |    NONE      |      NONE       |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`msaflagopen <msaflagopen_step>` |  MOS, IFU  |   MOS, IFU   |    MOS, IFU     |       MOS    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`extract_2d <extract_2d_step>`   |  MOS, FS   |   MOS, FS    |    MOS, FS      |       MOS    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`srctype <srctype_step>`         |  NONE      |    NONE      |      NONE       |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`wavecorr <wavecorr_step>`       |   ALL      |     ALL      |       ALL       |       ALL    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`flat_field <flatfield_step>`    |            |              |                 |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+
|              - D-FLAT                 |   ALL      |     ALL      |       ALL       |              |
+---------------------------------------+------------+--------------+-----------------+--------------+
|              - S-FLAT                 |   ALL      |    NONE      |       ALL       |              |
+---------------------------------------+------------+--------------+-----------------+--------------+
|              - F-FLAT                 |  NONE      |    NONE      |      NONE       |              |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`pathloss <pathloss_step>`       |  NONE      |    NONE      |      NONE       |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`barshadow <barshadow_step>`     |  NONE      |    NONE      |      NONE       |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`photom <photom_step>`           |  NONE      |    NONE      |      NONE       |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`resample_spec <resample_step>`  |  MOS, FS   |    NONE      |    MOS, FS      |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`cube_build <cube_build_step>`   |   IFU      |    NONE      |       IFU       |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+
| :ref:`extract_1d <extract_1d_step>`   |   ALL      |    NONE      |       ALL       |      NONE    |
+---------------------------------------+------------+--------------+-----------------+--------------+

In the :ref:`resample_spec <resample_step>` and :ref:`cube_build <cube_build_step>` steps, the spectra are
transformed to a space of (wavelength, offset along the slit) without applying a tangent plane projection.

Arguments
---------
The ``calwebb_spec2`` pipeline has two optional arguments.

``--save_bsub`` (boolean, default=False)
  If set to ``True``, the results of the background subtraction step will be saved
  to an intermediate file, using a product type of "_bsub" or "_bsubints", depending on
  whether the data are 2D (averaged over integrations) or 3D (per-integration results).

``--save_wfss_esec`` (boolean, default=False)
  If set to ``True``, an intermediate image product is created for WFSS exposures that
  is in units of electrons/sec, instead of the normal DN/sec units that are used throughout
  the rest of processing. This product can be useful for doing off-line specialized
  processing of WFSS images. This product is created after the :ref:`background <background_step>`
  and :ref:`flat-field <flatfield_step>` steps have been applied, but before the
  :ref:`extract_2d <extract_2d_step>` step, so that it is the full WFSS image. The conversion
  to units of electrons/sec is accomplished by loading the :ref:`GAIN <gain_reffile>` reference file,
  computing the mean gain across all pixels (excluding reference pixels), and multiplying the WFSS
  image by the mean gain.  The intermediate file will have a product type of "_esec".
  Only applies to WFSS exposures.

Inputs
------

2D or 3D countrate data
+++++++++++++++++++++++

:Data model: `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`,
             or `~jwst.datamodels.CubeModel`
:File suffix: _rate or _rateints

The input to the ``Spec2Pipeline`` pipeline is a countrate exposure, in the form
of either "_rate" or "_rateints" data. A single input file can be processed or an
ASN file listing multiple inputs can be used, in which case the processing steps
will be applied to each input exposure, one at a time.

If "_rateints" products are used as input, for modes other than NIRSpec Fixed Slit,
each step applies its algorithm to each integration in the exposure, where appropriate.
For the NIRSpec Fixed Slit mode the ``calwebb_spec2`` pipeline will currently
skip both the :ref:`resample_spec <resample_step>` step and the
:ref:`extract_1d <extract_1d_step>` step, because neither step supports
multiple integration input products for this mode.

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
+++++++++++++++++++++++++++++++++++

:Data model: `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`,
              or `~jwst.datamodels.CubeModel`
:File suffix: _bsub or _bsubints

This is an intermediate product that is only created if "--save_bsub" is set
to ``True`` and will contain the data as output from the :ref:`background <background_step>`
step. If the input is a "_rate" product, this will be a "_bsub" product, while
"_rateints" inputs will be saved as "_bsubints."

2D or 3D calibrated data
++++++++++++++++++++++++

:Data model: `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`,
             `~jwst.datamodels.CubeModel`,
             `~jwst.datamodels.SlitModel`, or `~jwst.datamodels.MultiSlitModel`
:File suffix: _cal or _calints

The output is a fully calibrated, but unrectified, exposure, using the product
type suffix "_cal" or "_calints", depending on the type of input,
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
the output type will consist of either a single slit model or a multi-slit model:

- NIRSpec Bright-Object and NIRCam TSO Grism: `~jwst.datamodels.SlitModel`
- NIRSpec Fixed Slit and MOS, as well as WFSS: `~jwst.datamodels.MultiSlitModel`

The multi-slit model is simply an array of multiple slit models, each one
containing the data and relevant meta data for a particular extracted slit or
source. A `~jwst.datamodels.MultiSlitModel` product will contain multiple
tuples of SCI, ERR, DQ, WAVELENGTH, etc. arrays; one for each of the
extracted slits/sources.

2D resampled data
+++++++++++++++++

:Data model: `~jwst.datamodels.SlitModel` or `~jwst.datamodels.MultiSlitModel`
:File suffix: _s2d

If the input is a 2D exposure type that gets resampled/rectified by the
:ref:`resample_spec <resample_step>` step, the rectified 2D spectral product is saved as a
"_s2d" file. This image is intended for use as a quick-look product only and is
not used in subsequent processing. The 2D unresampled, calibrated ("_cal")
products are passed along as input to subsequent Stage 3 processing.

If the input to the :ref:`resample_spec <resample_step>` step is a `~jwst.datamodels.MultiSlitModel`,
then the resampled output will be in the form of a
`~jwst.datamodels.MultiSlitModel`, which contains an array of individual models,
one per slit. Otherwise the output will be a single `~jwst.datamodels.SlitModel`.

3D resampled (IFU cube) data
++++++++++++++++++++++++++++

:Data model: `~jwst.datamodels.IFUCubeModel`
:File suffix: _s3d

If the data are NIRSpec IFU or MIRI MRS, the result of the :ref:`cube_build <cube_build_step>`
step will be 3D IFU spectroscopic cube saved to a "_s3d" file. The IFU cube is built from
the data contained in a single exposure and is intended for use as a quick-look
product only. The 2D unresampled, calibrated ("_cal") products are passed along as
input to subsequent Stage 3 processing.

1D extracted spectral data
++++++++++++++++++++++++++

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

.. include:: ../references_general/pars-spec2pipeline_reffile.inc
