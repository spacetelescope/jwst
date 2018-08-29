.. _pipelines:

Pipeline Modules
================

The pipelines that call individual correction steps in various
orders are defined as python classes within python code modules. The pipelines
can be executed by referencing their class name or through the use of a
configuration (.cfg) file that in turn references the class. The table below
shows the pipeline classes that are currently available, the
corresponding pre-defined configurations that make use of those classes, and
the instrument modes to which they can be applied.

+-------------------+-----------------------+------------------------------+
| Class Name        | Configuration File    | Used For                     |
+===================+=======================+==============================+
| Detector1Pipeline | calwebb_detector1.cfg | Stage 1: all non-TSO modes   |
+-------------------+-----------------------+------------------------------+
| Detector1Pipeline | calwebb_tso1.cfg      | Stage 1: all TSO modes       |
+-------------------+-----------------------+------------------------------+
| DarkPipeline      | calwebb_dark.cfg      | Stage 1: darks               |
+-------------------+-----------------------+------------------------------+
| GuiderPipeline    | calwebb_guider.cfg    | Stage 1+2: FGS guiding modes |
+-------------------+-----------------------+------------------------------+
| Image2Pipeline    | calwebb_image2.cfg    | Stage 2: imaging modes       |
+-------------------+-----------------------+------------------------------+
| Spec2Pipeline     | calwebb_spec2.cfg     | Stage 2: spectroscopy modes  |
+-------------------+-----------------------+------------------------------+
| Image3Pipeline    | calwebb_image3.cfg    | Stage 3: imaging modes       |
+-------------------+-----------------------+------------------------------+
| Spec3Pipeline     | calwebb_spec3.cfg     | Stage 3: spectroscopy modes  |
+-------------------+-----------------------+------------------------------+
| Ami3Pipeline      | calwebb_ami3.cfg      | Stage 3: NIRISS AMI mode     |
+-------------------+-----------------------+------------------------------+
| Coron3Pipeline    | calwebb_coron3.cfg    | Stage 3: Coronagraphic mode  |
+-------------------+-----------------------+------------------------------+
| TSO3Pipeline      | calwebb_tso3.cfg      | Stage 3: Time Series mode    |
+-------------------+-----------------------+------------------------------+

The data from different observing modes needs to be processed with
different combinations of the pipeline stages listed above. Observing
modes are usually identifiable via the value of the `EXP_TYPE` keyword in
the data product. The following table lists the pipeline modules that get
applied to each `EXP_TYPE` instance.

+---------------------+-------------------+------------------+------------------+
| | EXP_TYPE          | Stage 1 Pipeline  | Stage 2 Pipeline | Stage 3 Pipeline |
+=====================+===================+==================+==================+
| | FGS_IMAGE         | calwebb_detector1 | calwebb_image2   | calwebb_image3   |
+---------------------+-------------------+------------------+------------------+
| | FGS_FOCUS         | calwebb_detector1 | calwebb_image2   | N/A              |
+---------------------+-------------------+------------------+------------------+
| | FGS_DARK          | calwebb_dark1     | N/A              | N/A              |
+---------------------+-------------------+------------------+------------------+
| | FGS_SKYFLAT       | calwebb_detector1 | N/A              | N/A              |
| | FGS_INTFLAT       |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
|                     |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
| | MIR_IMAGE         | calwebb_detector1 | calwebb_image2   | calwebb_image3   |
+---------------------+-------------------+------------------+------------------+
| | MIR_MRS           | calwebb_detector1 | calwebb_spec2    | calwebb_spec3    |
+---------------------+-------------------+------------------+------------------+
| | MIR_LRS-FIXEDSLIT | calwebb_detector1 | calwebb_spec2    | calwebb_spec3    |
+---------------------+-------------------+------------------+------------------+
| | MIR_LRS-SLITLESS  | calwebb_tso1      | calwebb_spec2    | calwebb_tso3     |
+---------------------+-------------------+------------------+------------------+
| | MIR_LYOT          | calwebb_detector1 | calwebb_image2   | calwebb_coron3   |
| | MIR_4QPM          |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
| | MIR_TACQ          | calwebb_detector1 | calwebb_image2   | N/A              |
+---------------------+-------------------+------------------+------------------+
| | MIR_DARK          | calwebb_dark1     | N/A              | N/A              |
+---------------------+-------------------+------------------+------------------+
| | MIR_FLATIMAGE     | calwebb_detector1 | N/A              | N/A              |
| | MIR_FLATMRS       |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
|                     |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
| | NRC_IMAGE         | calwebb_detector1 | calwebb_image2   | calwebb_image3   |
+---------------------+-------------------+------------------+------------------+
| | NRC_CORON         | calwebb_detector1 | calwebb_image2   | calwebb_coron3   |
+---------------------+-------------------+------------------+------------------+
| | NRC_WFSS          | calwebb_detector1 | calwebb_spec2    | calwebb_spec3    |
+---------------------+-------------------+------------------+------------------+
| | NRC_TSIMAGE       | calwebb_tso1      | calwebb_image2   | calwebb_tso3     |
+---------------------+-------------------+------------------+------------------+
| | NRC_TSGRISM       | calwebb_tso1      | calwebb_spec2    | calwebb_tso3     |
+---------------------+-------------------+------------------+------------------+
| | NRC_TACQ          | calwebb_detector1 | calwebb_image2   | N/A              |
| | NRC_TACONFIRM     |                   |                  |                  |
| | NRC_FOCUS         |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
| | NRC_DARK          | calwebb_dark1     | N/A              | N/A              |
+---------------------+-------------------+------------------+------------------+
| | NRC_FLAT          | calwebb_detector1 | N/A              | N/A              |
| | NRC_LED           |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
|                     |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
| | NIS_IMAGE         | calwebb_detector1 | calwebb_image2   | calwebb_image3   |
+---------------------+-------------------+------------------+------------------+
| | NIS_WFSS          | calwebb_detector1 | calwebb_spec2    | calwebb_spec3    |
+---------------------+-------------------+------------------+------------------+
| | NIS_SOSS          | calwebb_tso1      | calwebb_spec2    | calwebb_tso3     |
+---------------------+-------------------+------------------+------------------+
| | NIS_AMI           | calwebb_detector1 | calwebb_image2   | calwebb_ami3     |
+---------------------+-------------------+------------------+------------------+
| | NIS_TACQ          | calwebb_detector1 | calwebb_image2   | N/A              |
| | NIS_TACONFIRM     |                   |                  |                  |
| | NIS_FOCUS         |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
| | NIS_DARK          | calwebb_dark1     | N/A              | N/A              |
+---------------------+-------------------+------------------+------------------+
| | NIS_LAMP          | calwebb_detector1 | N/A              | N/A              |
+---------------------+-------------------+------------------+------------------+
|                     |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
| | NRS_FIXEDSLIT     | calwebb_detector1 | calwebb_spec2    | calwebb_spec3    |
| | NRS_IFU           |                   |                  |                  |
| | NRS_MSASPEC       |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
| | NRS_BRIGHTOBJ     | calwebb_tso1      | calwebb_spec2    | calwebb_tso3     |
+---------------------+-------------------+------------------+------------------+
| | NRS_IMAGE         | calwebb_detector1 | calwebb_image2   | N/A              |
| | NRS_TACQ          |                   |                  |                  |
| | NRS_TACONFIRM     |                   |                  |                  |
| | NRS_BOTA          |                   |                  |                  |
| | NRS_TASLIT        |                   |                  |                  |
| | NRS_CONFIRM       |                   |                  |                  |
| | NRS_FOCUS         |                   |                  |                  |
| | NRS_MIMF          |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+
| | NRS_DARK          | calwebb_dark1     | N/A              | N/A              |
+---------------------+-------------------+------------------+------------------+
| | NRS_AUTOWAVE      | calwebb_detector1 | N/A              | N/A              |
| | NRS_AUTOFLAT      |                   |                  |                  |
| | NRS_LAMP          |                   |                  |                  |
+---------------------+-------------------+------------------+------------------+

Input Files, Output Files and Data Models
=========================================
An important concept used throughout the JWST pipeline is the :py:class:`Data
Model <jwst.datamodels.DataModel>`. Nearly all data used by any of the pipeline code is
encapsulated in a data model. Most input is read into a data model and
all output is produced by a data model. When possible, this document
will indicate the data model associated with a file type, usually as a
parenthetical link to the data model in question. For some steps, the
output file may represent different data models depending on the input
to those steps. As a result, the data models listed here will not be
an exhaustive list.

.. _stage1-flow:

Stage 1 Pipeline Step Flow (calwebb_detector1)
==============================================
Stage 1 processing applies basic detector-level corrections to all exposure
types (imaging, spectroscopic, coronagraphic, etc.). It is applied to one
exposure at a time. The pipeline module for stage 1 processing is
``calwebb_detector1`` (the equivalent pipeline class is ``Detector1Pipeline``). It is
often referred to as ``ramps-to-slopes`` processing, because the input raw data
are in the form of one or more ramps (integrations) containing accumulating
counts from the non-destructive detector readouts and the output is a corrected
countrate (slope) image. The list of steps applied by the Build 7.1 calwebb_detector1
pipeline is as follows.

================= =================
calwebb_detector1 calwebb_detector1
(All Near-IR)     (MIRI)
================= =================
group_scale       group_scale
dq_init           dq_init
saturation        saturation
ipc               ipc
superbias         linearity
refpix            rscd
linearity         lastframe
persistence       dark_current
dark_current      refpix
\                 persistence
jump              jump
ramp_fit          ramp_fit
gain_scale        gain_scale
================= =================

If the ``calwebb_tso1.cfg`` configuration file is used to execute this pipeline,
the ``ipc``, ``lastframe``, and ``persistence`` steps will be skipped.

Inputs
------

* Raw 4D product: The input to ``calwebb_detector1`` is a single raw exposure file,
  e.g. ``jw80600012001_02101_00003_mirimage_uncal.fits``, which contains the
  original raw data from all of the detector readouts in the exposure
  (ncols x nrows x ngroups x nintegrations).

Outputs
-------

* 2D Countrate product: All types of inputs result in a 2D countrate product,
  resulting from averaging over all of the integrations within the exposure.
  The output file will be of type ``_rate``, e.g.
  ``jw80600012001_02101_00003_mirimage_rate.fits``.

* 3D Countrate product: If the input exposure contains more than one integration
  (NINTS>1), a 3D countrate product is created that contains the individual
  results of each integration. The 2D countrate images for each integration are
  stacked along the 3rd axis of the data cubes (ncols x nrows x nints). This
  output file will be of type ``_rateints``.

Arguments
---------
The ``calwebb_detector1`` pipeline has one optional argument:

* ``save_calibrated_ramp``

which is a boolean argument with a default value of ``False``. If the user sets
it to ``True``, the pipeline will save intermediate data to a file as it
exists at the end of the ``jump`` step (just before ramp fitting). The data at
this stage of the pipeline are still in the form of the original 4D ramps
(ncols x nrows x ngroups x nints) and have had all of the detector-level
correction steps applied to it, including the detection and flagging of
Cosmic-Ray hits within each ramp (integration). If created, the name of the
intermediate file will be constructed from the root name of the input file, with
the new product type suffix ``_ramp`` appended
(e.g. ``jw80600012001_02101_00003_mirimage_ramp.fits``).

Dark Pipeline Step Flow (calwebb_dark)
======================================
The stage 1 dark (``calwebb_dark``) processing pipeline is intended for use
with dark exposures. It applies all of the same detector-level correction steps
as the ``calwebb_detector1`` pipeline, but stops just before the application of the
``dark_current`` step.

Inputs
------

* Raw 4D Dark product: The input to ``calwebb_dark`` is a single raw dark
  exposure.

Outputs
-------

* 4D Corrected product: The output is a 4D (ncols x nrows x ngroups x nints)
  product that has had all corrections up to, but not including, the
  ``dark_current`` step, with a product file type of ``_dark``.

Arguments
---------
The ``calwebb_dark`` pipeline does not have any optional arguments.

Guider Pipeline Step Flow (calwebb_guider)
==========================================
The guider (``calwebb_guider``) processing pipeline is only for use with FGS
guiding mode exposures (ID, ACQ1, ACQ2, TRACK, and FineGuide).
It applies three detector-level correction and calibration steps to uncalibrated
guider data files, as listed in the table below.

+----------------+
| calwebb_guider |
+================+
| dq_init        |
+----------------+
| guider_cds     |
+----------------+
| flat_field     |
+----------------+

Inputs
------

* Raw 4D guide-mode product: The input to ``calwebb_guider`` is a single raw
  guide-mode data file.

Outputs
-------

* 3D Calibrated product: The output is a 3D (ncols x nrows x nints)
  countrate product that has been flat-fielded and has bad pixels flagged.
  See the documentation for the guider_cds step for details on the
  conversion from raw readouts to countrate images.

Arguments
---------
The ``calwebb_guider`` pipeline does not have any optional arguments.

.. _stage2-imaging-flow:

Stage 2 Imaging Pipeline Step Flow (calwebb_image2)
====================================================
Stage 2 imaging (``calwebb_image2``) processing applies additonal corrections
that result in a fully calibrated individual exposure. The list of correction
steps applied by the calwebb_image2 imaging pipeline is as follows.

+----------------+
| calwebb_image2 |
+================+
| background     |
+----------------+
| assign_wcs     |
+----------------+
| flat_field     |
+----------------+
| photom         |
+----------------+
| resample       |
+----------------+

Inputs
------

* 2D or 3D Countrate product: The input to the ``calwebb_image2`` pipeline is
  a countrate exposure, in the form of either ``_rate`` or ``_rateints``
  files. A single input file can be processed or an ASN file listing
  multiple inputs can be used, in which case the processing steps will be
  applied to each input exposure, one at a time. If ``_rateints`` products are
  used as input, the steps in the pipeline are applied individually to each
  integration in an exposure, where appropriate.

Outputs
-------

* 2D or 3D Calibrated product: The output is a calibrated exposure, using
  the product type suffix ``_cal`` or ``_calints``, depending on the type of
  input (e.g. ``jw80600012001_02101_00003_mirimage_cal.fits``).

Arguments
---------
The ``calwebb_image2`` pipeline has one optional argument ``save_bsub``,
which is set to ``False`` by default. If set to ``True``, the results of
the background subtraction step will be saved to an intermediate file,
using a product type of ``_bsub`` or ``_bsubints`` (depending on the type
of input).

.. _stage2-spectroscopic-flow:

Stage 2 Spectroscopic Pipeline Step Flow (calwebb_spec2)
==========================================================
Stage 2 spectroscopic (``calwebb_spec2``) pipeline applies additional
corrections to countrate products that result in fully calibrated individual
exposures.
The list of correction steps is shown below. Some steps are only applied to
certain instruments or instrument modes, as noted in the table.

+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| Instrument Mode      |     NIRSpec    |     MIRI      |    NIRISS   | NIRCam |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| Step                 | FS | MOS | IFU | FS | SL | MRS | SOSS | WFSS | WFSS   |
+======================+====+=====+=====+====+====+=====+======+======+========+
| assign_wcs           | X  |  X  |  X  | X  | X  |  X  |   X  |  X   |   X    |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| background           | X  |  X  |  X  | X  | X  |  X  |   X  |  X   |   X    |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| imprint              |    |  X  |  X  |    |    |     |      |      |        |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| msaflagopen          |    |  X  |  X  |    |    |     |      |      |        |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| extract_2d\ :sup:`1` | X  |  X  |     |    |    |     |      |  X   |   X    |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| flat_field\ :sup:`1` | X  |  X  |  X  | X  | X  |  X  |   X  |  X   |   X    |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| srctype              | X  |  X  |  X  | X  | X  |  X  |   X  |  X   |   X    |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| straylight           |    |     |     |    |    |  X  |      |      |        |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| fringe               |    |     |     |    |    |  X  |      |      |        |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| pathloss             | X  |  X  |  X  |    |    |     |      |      |        |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| barshadow            |    |  X  |     |    |    |     |      |      |        |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| photom               | X  |  X  |  X  | X  | X  |  X  |   X  |  X   |   X    |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| resample_spec        | X  |  X  |     |    |    |     |      |      |        |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| cube_build           |    |     |  X  |    |    |  X  |      |      |        |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+
| extract_1d           | X  |  X  |  X  | X  | X  |  X  |   X  |  X   |   X    |
+----------------------+----+-----+-----+----+----+-----+------+------+--------+

:sup:`1`\ Note that the order of the extract_2d and flat_field steps is reversed
(flat_field is performed first) for NIRISS and NIRCam WFSS exposures.

The ``resample_spec`` step produces a resampled/rectified product for non-IFU
modes of some spectroscopic exposures. If the ``resample_spec`` step
is not applied to a given exposure, the ``extract_1d`` operation will be
performed on the original (unresampled) data. The ``cube_build`` step produces
a resampled/rectified cube for IFU exposures, which is then used as input to
the ``extract_1d`` step.

Inputs
------
The input to the ``calwebb_spec2`` pipeline can be either a single countrate
(``_rate`` or ``_rateints``) exposure or an Association (ASN) file
listing multiple exposures. The background subtraction (``bkg_subtract``) and
imprint subtraction (``imprint_subtract``) steps can only be executed when
the pipeline is supplied with an association of exposures, because they rely
on multiple exposures to perform their tasks. The ASN file must not only list
the input exposures, but must also contain information that indicates their
relationships to one another.

The background subtraction step can be applied to an assocation containing
nodded exposures, such as for MIRI LRS fixed-slit, NIRSpec fixed-slit, and
NIRSpec MSA observations, or an association that contains dedicated exposures
of a background target. The step will accomplish background subtraction by
doing direct subtraction of nodded exposures from one another or by direct
subtraction of dedicated background expsoures from the science target exposures.

Background subtraction for Wide-Field Slitless Spectroscopy (WFSS) exposures
is accomplished by scaling and subtracting a master background image from a
CRDS reference file.

The imprint subtraction step, which is only applied to NIRSpec MSA and IFU
exposures, also requires the use of an ASN file, in order to specify which of
the inputs is to be used as the imprint exposure. The imprint exposure will be
subtracted from all other exposures in the association.

If a single countrate product is used as input, the background subtraction
and imprint subtraction steps will be skipped and only the remaining regular
calibration steps will be applied to the input exposure.

Outputs
-------
Two or three different types of outputs are created by ``calwebb_spec2``.

* Calibrated product: All types of inputs result in a fully-calibrated
  product at the end of the ``photom`` step, which uses the ``_cal`` or
  ``_calints`` product type suffix, depending on whether the input was a
  ``_rate`` or ``_rateints`` product, respectively.

* Resampled 2D product: If the input is a 2D exposure type that gets
  resampled/rectified by the ``resample_spec`` step, the rectified 2D spectral
  product created by the ``resample_spec`` step is saved as a ``_s2d`` file.
  3D (``_rateints``) input exposures do not get resampled.

* Resampled 3D product: If the data are NIRSpec IFU or MIRI MRS, the
  result of the ``cube_build`` step will be saved as a ``_s3d`` file.

* 1D Extracted Spectrum product: All types of inputs result in a 1D extracted
  spectral data product, which is saved as a ``_x1d`` or ``_x1dints`` file,
  depending on the input type.

If the input to ``calwebb_spec2`` is an ASN file, these products are created
for each input exposure.

Arguments
---------
The ``calwebb_spec2`` pipeline has one optional argument:

* ``save_bsub``

which is a Boolean argument with a default value of ``False``. If the user sets
it to ``True``, the results of the background subtraction step (if applied) are
saved to an intermediate file of type ``_bsub`` or ``_bsubints``, as appropriate.

.. _stage3-imaging-flow:

Stage 3 Imaging Pipeline Step Flow (calwebb_image3)
===================================================
Stage 3 processing for imaging observations is intended for combining the 
calibrated data
from multiple exposures (e.g. a dither or mosaic pattern) into a single
rectified (distortion corrected) product.
Before being combined, the exposures receive additional corrections for the
purpose of astrometric alignment, background matching, and outlier rejection.
The steps applied by the ``calwebb_image3`` pipeline are shown below.

+-------------------+
| calwebb_image3    |
+===================+
| tweakreg          |
+-------------------+
| skymatch          |
+-------------------+
| outlier_detection |
+-------------------+
| resample          |
+-------------------+
| source_catalog    |
+-------------------+

Inputs
------

* Associated 2D Calibrated products: The inputs to ``calwebb_image3`` will
  usually be in the form of an ASN file that lists multiple exposures to be
  processed and combined into a single output product. The individual exposures
  should be calibrated (``_cal``) products from ``calwebb_image2`` processing.

* Single 2D Calibrated product: It is also possible use a single ``_cal`` file
  as input to ``calwebb_image3``, in which case only the ``resample`` and
  ``source_catalog`` steps will be applied.

Outputs
-------

* Resampled 2D Image product (:py:class:`DrizProductModel
  <jwst.datamodels.DrizProductModel>`): A resampled/rectified 2D image product of type
  ``_i2d`` is created containing the rectified single exposure or the rectified
  and combined association of exposures, which is the direct output of the
  ``resample`` step.

* Source catalog: A source catalog produced from the ``_i2d`` product is saved
  as an ASCII file in ``ecsv`` format, with a product type of ``_cat``.

* CR-flagged products: If the ``outlier_detection`` step is applied, a new version
  of each input calibrated exposure product is created, which contains a DQ array
  that has been updated to flag pixels detected as outliers. This updated
  product is known as a CR-flagged product and the file is identified by including
  the association candidate ID in the original input ``_cal`` file name and
  changing the product type to ``_crf``, e.g.
  ``jw96090001001_03101_00001_nrca2_o001_crf.fits``.

.. _stage3-spectroscopic-flow:

Stage 3 Spectroscopic Pipeline Step Flow (calwebb_spec3)
=========================================================
Stage 3 processing for spectroscopic observations is intended for combining the 
calibrated data from multiple exposures (e.g. a dither pattern) into a single
rectified (distortion corrected) product and a combined 1D spectrum.
Before being combined, the exposures may receive additional corrections for the
purpose of background matching and outlier rejection.
The steps applied by the ``calwebb_spec3`` pipeline are shown below.

+-------------------+----+-----+-----+----+-----+--------+--------+
| Instrument Mode   |     NIRSpec    |   MIRI   | NIRISS | NIRCam |
+-------------------+----+-----+-----+----+-----+--------+--------+
| Step              | FS | MOS | IFU | FS | MRS | WFSS   | WFSS   |
+===================+====+=====+=====+====+=====+========+========+
| mrs_imatch        |    |     |     |    |  X  |        |        |
+-------------------+----+-----+-----+----+-----+--------+--------+
| outlier_detection | X  |  X  |  X  | X  |  X  |   X    |   X    |
+-------------------+----+-----+-----+----+-----+--------+--------+
| resample_spec     | X  |  X  |     | X  |     |   X    |   X    |
+-------------------+----+-----+-----+----+-----+--------+--------+
| cube_build        |    |     |  X  |    |  X  |        |        |
+-------------------+----+-----+-----+----+-----+--------+--------+
| extract_1d        | X  |  X  |  X  | X  |  X  |   X    |   X    |
+-------------------+----+-----+-----+----+-----+--------+--------+

NOTE: In B7.1 the ``calwebb_spec3`` pipeline is very much a prototype and
not all steps are functioning properly for all modes. In particular, the
``outlier_detection`` step does not yet work well, if at all, for any of
the spectroscopic modes. Also, the ``resample_spec`` step does not work
for dithered slit-like spectra (i.e. all non-IFU modes). Processing of
NIRSpec IFU and MIRI MRS exposures does work, using the
``mrs_imatch``, ``cube_build``, and ``extract_1d`` steps.

Inputs
------

* Associated 2D Calibrated products: The inputs to ``calwebb_spec3`` will
  usually be in the form of an ASN file that lists multiple exposures to be
  processed and combined into a single output product. The individual exposures
  should be calibrated (``_cal``) products from ``calwebb_spec2`` processing.

Outputs
-------

* CR-flagged products: If the ``outlier_detection`` step is applied, a new version
  of each input calibrated exposure product is created, which contains a DQ array
  that has been updated to flag pixels detected as outliers. This updated
  product is known as a CR-flagged product and the file is identified by including
  the association candidate ID in the original input ``_cal`` file name and
  changing the product type to ``_crf``, e.g.
  ``jw96090001001_03101_00001_nrs2_o001_crf.fits``.

* Resampled 2D spectral product (:py:class:`DrizProductModel
  <jwst.datamodels.DrizProductModel>`): A resampled/rectified 2D product of type
  ``_s2d`` is created containing the rectified and combined association of
  exposures, which is the direct output of the ``resample_spec`` step, when
  processing non-IFU modes.

* Resampled 3D spectral product (:py:class:`IFUCubeModel
  <jwst.datamodels.IFUCubeModel>`): A resampled/rectified 3D product of type
  ``_s3d`` is created containing the rectified and combined association of
  exposures, which is the direct output of the ``cube_build`` step, when
  processing IFU modes.

* 1D Extracted Spectrum product: All types of inputs result in a 1D extracted
  spectral data product, which is saved as a ``_x1d`` file, and is the result
  of performing 1D extraction on the combined ``_s2d`` or ``_s3d`` product.

.. _stage3-ami-flow:

Stage 3 Aperture Masking Interferometry (AMI) Pipeline Step Flow (calwebb_ami3)
===============================================================================
The stage 3 AMI (``calwebb_ami3``) pipeline is to be applied to
associations of calibrated NIRISS AMI exposures and is used to compute fringe
parameters and correct science target fringe parameters using observations of
reference targets.
The steps applied by the ``calwebb_ami3`` pipeline are shown below.

+---------------+
| calwebb_ami3  |
+===============+
| ami_analyze   |
+---------------+
| ami_average   |
+---------------+
| ami_normalize |
+---------------+

Inputs
------

* Associated 2D Calibrated products: The inputs to ``calwebb_ami3`` need to be
  in the form of an ASN file that lists multiple science target exposures,
  and optionally reference target exposures as well. The individual exposures
  should be in the form of calibrated (``_cal``) products from ``calwebb_image2``
  processing.

Outputs
-------

* AMI product (:py:class:`AmiLgModel <jwst.datamodels.AmiLgModel>`):
  For every input exposure, the fringe parameters and closure phases caculated
  by the ``ami_analyze`` step are saved to an ``_ami`` product file, which is
  a table containing the fringe parameters and closure phases. Product names
  use the original input ``_cal`` file name, with the association candidate ID
  included and the product type changed to ``_ami``, e.g.
  ``jw93210001001_03101_00001_nis_a0003_ami.fits``.

* Averaged AMI product (:py:class:`AmiLgModel <jwst.datamodels.AmiLgModel>`):
  The AMI results averaged over all science or reference
  exposures, calculated by the ``ami_average`` step, are saved to an ``_amiavg``
  product file. Separate products are created for the science target and
  reference target data. Note that these output products are only created if the
  pipeline argument ``save_averages`` (see below) is set to ``True``.

* Normalized AMI product (:py:class:`AmiLgModel <jwst.datamodels.AmiLgModel>`):
  If reference target exposures are included in the input
  ASN, the averaged AMI results for the science target will be normalized by the
  averaged AMI results for the reference target, via the ``ami_normalize`` step,
  and will be saved to an ``_aminorm`` product file.

Arguments
---------
The ``calwebb_ami3`` pipeline has one optional argument:

* ``save_averages``

which is a Boolean parameter set to a default value of ``False``. If the user
sets this agument to ``True``, the results of the ``ami_average`` step will be
saved, as described above.

.. _stage3-coron-flow:

Stage 3 Coronagraphic Pipeline Step Flow (calwebb_coron3)
===============================================================================
The stage 3 coronagraphic (``calwebb_coron3``) pipeline is to be applied to
associations of calibrated NIRCam coronagraphic and MIRI Lyot and 4QPM
exposures, and is used to produce psf-subtracted, resampled, combined images
of the source object.

The steps applied by the ``calwebb_coron3`` pipeline are shown in the table
below.

+----------------------------------------------------------------------------------------------------+
| :py:class:`calwebb_coron3 <jwst.pipeline.calwebb_coron3.Coron3Pipeline>`                           |
+====================================================================================================+
| :py:class:`stack_refs <jwst.coron.stack_refs_step.StackRefsStep>`                                  |
+----------------------------------------------------------------------------------------------------+
| :py:class:`align_refs <jwst.coron.align_refs_step.AlignRefsStep>`                                  |
+----------------------------------------------------------------------------------------------------+
| :py:class:`klip <jwst.coron.klip_step.KlipStep>`                                                   |
+----------------------------------------------------------------------------------------------------+
| :py:class:`outlier_detection <jwst.outlier_detection.outlier_detection_step.OutlierDetectionStep>` |
+----------------------------------------------------------------------------------------------------+
| :py:class:`resample <jwst.resample.resample_step.ResampleStep>`                                    |
+----------------------------------------------------------------------------------------------------+


Inputs
------

* Associated Calibrated products: The input to ``calwebb_coron3`` is assumed
  to be in the form of an ASN file that lists multiple observations of
  a science target and, optionally, a reference PSF target. The individual science
  target and PSF reference exposures should be in the form of 3D calibrated (``_calints``)
  products from ``calwebb_image2`` processing.

Outputs
-------

* Stacked PSF images: The data from each input PSF reference exposure are
  concatenated into a single combined 3D stack, for use by subsequent steps. The
  stacked PSF data gets written to disk in the form of a psfstack (``_psfstack``)
  product from
  :py:class:`stack_refs step <jwst.coron.stack_refs_step.StackRefsStep>`.

* Aligned PSF images: The initial processing requires aligning all input PSFs
  specified in the ASN.  The aligned PSF images then gets written to disk in the
  form of psfalign (``_psfalign``) products from
  :py:class:`align_refs step <jwst.coron.align_refs_step.AlignRefsStep>`.

* PSF-subtracted exposures: The :py:class:`klip step <jwst.coron.klip_step.KlipStep>`
  takes the aligned PSF images and applies them to each of the science exposures
  in the ASN to create psfsub (``_psfsub``) products.

* CR-flagged products: The
  :py:class:`~jwst.outlier_detection.outlier_detection_step.OutlierDetectionStep`
  step is applied to the psfsub products to flag pixels in the DQ array
  that have been detected as outliers. This updated
  product is known as a CR-flagged product. A outlier-flagged product of
  type ``_crfints`` is created and can optionally get written to disk.

* Resampled product: The
  :py:class:`resample step <jwst.resample.resample_step.ResampleStep>` is
  applied to the CR-flagged products to create a single resampled, combined
  product for the science target.  This resampled product of type ``_i2d`` gets
  written to disk and returned as the final product from this pipeline.

.. _stage3-tso-flow:

Stage 3 Time-Series Observation(TSO) Pipeline Step Flow (calwebb_tso3)
===============================================================================
The stage 3 TSO (``calwebb_tso3``) pipeline is to be applied to
associations of calibrated TSO exposures (NIRCam TS imaging, NIRCam TS grism,
NIRISS SOSS, NIRSpec BrightObj, MIRI LRS Slitless) and is used to
produce calibrated time-series photometry of the source object.

The steps applied by the ``calwebb_tso3`` pipeline for an Imaging TSO observation
are shown below:

+----------------------------------------------------------------------------------------------------+
| :py:class:`calwebb_tso3 <jwst.pipeline.calwebb_tso3.Tso3Pipeline>`                                 |
+====================================================================================================+
| :py:class:`outlier_detection <jwst.outlier_detection.outlier_detection_step.OutlierDetectionStep>` |
+----------------------------------------------------------------------------------------------------+
| :py:class:`tso_photometry <jwst.tso_photometry.tso_photometry_step.TSOPhotometryStep>`             |
+----------------------------------------------------------------------------------------------------+

The steps applied by the ``calwebb_tso3`` pipeline for a Spectroscopic TSO
observation are shown below:

+----------------------------------------------------------------------------------------------------+
| :py:class:`calwebb_tso3 <jwst.pipeline.calwebb_tso3.Tso3Pipeline>`                                 |
+====================================================================================================+
| :py:class:`outlier_detection <jwst.outlier_detection.outlier_detection_step.OutlierDetectionStep>` |
+----------------------------------------------------------------------------------------------------+
| :py:class:`extract_1d <jwst.extract_1d.extract_1d_step.Extract1dStep>`                             |
+----------------------------------------------------------------------------------------------------+
| :py:class:`white_light <jwst.white_light.white_light_step.WhiteLightStep>`                         |
+----------------------------------------------------------------------------------------------------+

Inputs
------

* Associated 3D Calibrated products: The input to ``calwebb_tso3`` is assumed
  to be in the form of an ASN file that lists multiple science observations of
  a science target. The individual exposures should be in the form of 3D calibrated
  (``_calints``) products from either ``calwebb_image2`` or ``calwebb_spec2``
  processing.

Outputs
-------

* CR-flagged products: If the
  :py:class:`~outlier_detection.outlier_detection_step.OutlierDetectionStep`
  step is applied, a new version
  of each input calibrated exposure product is created, which contains a DQ array
  that has been updated to flag pixels detected as outliers. This update
  product is known as a CR-flagged product. A outlier-flagged product of
  type ``_crfints`` is created and can optionally get written to disk.

* Source photometry catalog for imaging TS observations: A source catalog produced
  from the ``_crfints`` product is saved as an ASCII file in ``ecsv`` format
  with a product type of ``_phot``.

* Extracted 1D spectra for spectroscopic TS observations: The ``extract_1d`` step is
  applied to create a ``MultiSpecModel`` for the entire set of
  spectra, with a product type of ``_x1dints``.

* White-light photometry for spectroscopic TS observations:  The ``white_light`` step
  is applied to the ``_x1dints`` extracted data to produce an ASCII catalog
  in ``ecsv`` format with a product type of ``_whtlt``, containing
  the wavelength-integrated white-light photometry of the source object.
