.. _pipelines:

Pipeline Stages
===============

End-to-end calibration of JWST data is divided into 3 main stages of
processing:

- Stage 1 consists of detector-level corrections that are performed on a
  group-by-group basis, followed by ramp fitting. The output of stage 1
  processing is a countrate image per exposure, or per integration for
  some modes. Details of this pipeline can be found at:

  - :ref:`calwebb_detector1`

- Stage 2 processing consists of additional instrument-level and
  observing-mode corrections and calibrations to produce fully calibrated
  exposures. The details differ for imaging and spectroscopic exposures,
  and there are some corrections that are unique to certain instruments or modes.
  Details are at:

  - :ref:`calwebb_image2`
  - :ref:`calwebb_spec2`

- Stage 3 processing consists of routines that work with multiple exposures
  and in most cases produce some kind of combined product.
  There are unique pipeline modules for stage 3 processing of
  imaging, spectroscopic, coronagraphic, AMI, and TSO observations. Details
  of each are available at:

  - :ref:`calwebb_image3`
  - :ref:`calwebb_spec3`
  - :ref:`calwebb_coron3`
  - :ref:`calwebb_ami3`
  - :ref:`calwebb_tso3`

In addition, there are several pipeline modules designed for special instrument or
observing modes, including:

- :ref:`calwebb_dark <calwebb_dark>` for processing dark exposures
- :ref:`calwebb_guider <calwebb_guider>` for calibrating FGS guide star data
- :ref:`calwebb_wfs-image2 <calwebb_wfs-image2>` for stage 2 WFS&C processing
- :ref:`calwebb_wfs-image3 <calwebb_wfs-image3>` for stage 3 WFS&C processing

Pipeline Classes and Configuration Files
========================================

Each pipeline consists of a certain sequence of calibration steps and is
defined as a python class within a python code module. The pipelines
can be executed from the command line either by referencing their class name or
by supplying a configuration (.cfg) file that in turn references the pipeline class
(see :ref:`strun_command_line`).
From within python, the pipelines are called by their class names, but
configuration files can still be supplied in order to set pipeline or step
parameter values (see :ref:`run_from_python`).
The table below lists the pipeline classes that are currently available, the
corresponding configuration files that call those classes, and
the observing modes for which they are intended. Note that there are some
pipeline modules that are referenced by more than one configuration file.

+------------------------------------+---------------------------+------------------------------+
| Pipeline Class                     | Configuration File        | Used For                     |
+====================================+===========================+==============================+
| `~jwst.pipeline.Detector1Pipeline` | calwebb_detector1.cfg     | Stage 1: all non-TSO modes   |
+                                    +---------------------------+------------------------------+
|                                    | calwebb_tso1.cfg          | Stage 1: all TSO modes       |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.pipeline.DarkPipeline`      | calwebb_dark.cfg          | Stage 1: darks               |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.pipeline.GuiderPipeline`    | calwebb_guider.cfg        | Stage 1+2: FGS guiding modes |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.pipeline.Image2Pipeline`    | calwebb_image2.cfg        | Stage 2: imaging modes       |
+                                    +---------------------------+------------------------------+
|                                    | calwebb_tso-image2.cfg    | Stage 2: TSO imaging modes   |
+                                    +---------------------------+------------------------------+
|                                    | calwebb_wfs-image2.cfg    | Stage 2: WFS&C imaging       |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.pipeline.Spec2Pipeline`     | calwebb_spec2.cfg         | Stage 2: spectroscopy modes  |
+                                    +---------------------------+------------------------------+
|                                    | calwebb_tso-spec2.cfg     | Stage 2: TSO spectral modes  |
+                                    +---------------------------+------------------------------+
|                                    | calwebb_nrslamp-spec2.cfg | Stage 2: NIRSpec lamps       |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.pipeline.Image3Pipeline`    | calwebb_image3.cfg        | Stage 3: imaging modes       |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.wfs_combine.WfsCombineStep` | calwebb_wfs-image3.cfg    | Stage 3: WFS&C imaging       |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.pipeline.Spec3Pipeline`     | calwebb_spec3.cfg         | Stage 3: spectroscopy modes  |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.pipeline.Ami3Pipeline`      | calwebb_ami3.cfg          | Stage 3: NIRISS AMI mode     |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.pipeline.Coron3Pipeline`    | calwebb_coron3.cfg        | Stage 3: Coronagraphic mode  |
+------------------------------------+---------------------------+------------------------------+
| `~jwst.pipeline.Tso3Pipeline`      | calwebb_tso3.cfg          | Stage 3: TSO modes           |
+------------------------------------+---------------------------+------------------------------+

Pipelines vs. Exposure Type
===========================

The data from different observing modes needs to be processed with
different combinations of the pipeline stages listed above. The proper pipeline
selection is usually based solely on the exposure type (EXP_TYPE keyword value).
Some modes, however, require additional selection criteria, such as whether the
data are to be treated as Time-Series Observations (TSO). Some EXP_TYPEs are
exclusively TSO, while others depend on the value of the TSOVISIT keyword.
The following table lists the pipeline modules that should get applied to various
observing modes, based on these selectors. Exposure types that do not allow TSO
mode are marked as "N/A" in the TSOVISIT column.

+---------------------+----------+-------------------+-----------------------+------------------+
| | EXP_TYPE          | TSOVISIT | Stage 1 Pipeline  | Stage 2 Pipeline      | Stage 3 Pipeline |
+=====================+==========+===================+=======================+==================+
| | FGS_DARK          | N/A      | calwebb_dark      | N/A                   | N/A              |
+---------------------+----------+-------------------+-----------------------+------------------+
| | FGS_SKYFLAT       | N/A      | calwebb_detector1 | N/A                   | N/A              |
| | FGS_INTFLAT       |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | FGS_FOCUS         | N/A      | calwebb_detector1 | calwebb_image2        | N/A              |
+---------------------+----------+-------------------+-----------------------+------------------+
| | FGS_IMAGE         | N/A      | calwebb_detector1 | calwebb_image2        | calwebb_image3   |
+---------------------+----------+-------------------+-----------------------+------------------+
| | FGS_ID-STACK      | N/A      | calwebb_guider    | N/A                   | N/A              |
| | FGS_ID-IMAGE      |          |                   |                       |                  |
| | FGS_ACQ1          |          |                   |                       |                  |
| | FGS_ACQ2          |          |                   |                       |                  |
| | FGS_TRACK         |          |                   |                       |                  |
| | FGS_FINEGUIDE     |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
|                     |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | MIR_DARKIMG       | N/A      | calwebb_dark      | N/A                   | N/A              |
| | MIR_DARKMRS       |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | MIR_FLATIMAGE     | N/A      | calwebb_detector1 | N/A                   | N/A              |
| | MIR_FLATIMAGE-EXT |          |                   |                       |                  |
| | MIR_FLATMRS       |          |                   |                       |                  |
| | MIR_FLATMRS-EXT   |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | MIR_TACQ          | N/A      | calwebb_detector1 | calwebb_image2        | N/A              |
+---------------------+----------+-------------------+-----------------------+------------------+
| | MIR_CORONCAL      | N/A      | calwebb_detector1 | calwebb_image2        | N/A              |
+---------------------+----------+-------------------+-----------------------+------------------+
| | MIR_IMAGE         | False    | calwebb_detector1 | calwebb_image2        | calwebb_image3   |
+                     +----------+-------------------+-----------------------+------------------+
|                     | True     | calwebb_tso1      | calwebb_tso-image2    | calwebb_tso3     |
+---------------------+----------+-------------------+-----------------------+------------------+
| | MIR_LRS-FIXEDSLIT | N/A      | calwebb_detector1 | calwebb_spec2         | calwebb_spec3    |
+---------------------+----------+-------------------+-----------------------+------------------+
| | MIR_LRS-SLITLESS  | True     | calwebb_tso1      | calwebb_tso-spec2     | calwebb_tso3     |
+                     +----------+-------------------+-----------------------+------------------+
|                     | False    | calwebb_detector1 | calwebb_spec2         | N/A              |
+---------------------+----------+-------------------+-----------------------+------------------+
| | MIR_MRS           | N/A      | calwebb_detector1 | calwebb_spec2         | calwebb_spec3    |
+---------------------+----------+-------------------+-----------------------+------------------+
| | MIR_LYOT          | N/A      | calwebb_detector1 | calwebb_image2        | calwebb_coron3   |
| | MIR_4QPM          |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
|                     |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRC_DARK          | N/A      | calwebb_dark      | N/A                   | N/A              |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRC_FLAT          | N/A      | calwebb_detector1 | N/A                   | N/A              |
| | NRC_LED           |          |                   |                       |                  |
| | NRC_GRISM         |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRC_TACQ          | N/A      | calwebb_detector1 | calwebb_image2        | N/A              |
| | NRC_TACONFIRM     |          |                   |                       |                  |
| | NRC_FOCUS         |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRC_IMAGE         | N/A      | calwebb_detector1 | calwebb_image2        | calwebb_image3   |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRC_CORON         | N/A      | calwebb_detector1 | calwebb_image2        | calwebb_coron3   |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRC_WFSS          | N/A      | calwebb_detector1 | calwebb_spec2         | calwebb_spec3    |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRC_TSIMAGE       | True     | calwebb_tso1      | calwebb_tso-image2    | calwebb_tso3     |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRC_TSGRISM       | True     | calwebb_tso1      | calwebb_tso-spec2     | calwebb_tso3     |
+---------------------+----------+-------------------+-----------------------+------------------+
|                     |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NIS_DARK          | N/A      | calwebb_dark      | N/A                   | N/A              |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NIS_LAMP          | N/A      | calwebb_detector1 | N/A                   | N/A              |
| | NIS_EXTCAL        |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NIS_TACQ          | N/A      | calwebb_detector1 | calwebb_image2        | N/A              |
| | NIS_TACONFIRM     |          |                   |                       |                  |
| | NIS_FOCUS         |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NIS_IMAGE         | N/A      | calwebb_detector1 | calwebb_image2        | calwebb_image3   |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NIS_AMI           | N/A      | calwebb_detector1 | calwebb_image2        | calwebb_ami3     |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NIS_WFSS          | N/A      | calwebb_detector1 | calwebb_spec2         | calwebb_spec3    |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NIS_SOSS          | True     | calwebb_tso1      | calwebb_tso-spec2     | calwebb_tso3     |
+                     +----------+-------------------+-----------------------+------------------+
|                     | False    | calwebb_detector1 | calwebb_spec2         | calwebb_spec3    |
+---------------------+----------+-------------------+-----------------------+------------------+
|                     |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRS_DARK          | N/A      | calwebb_dark      | N/A                   | N/A              |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRS_AUTOWAVE      | N/A      | calwebb_detector1 | calwebb_nrslamp-spec2 | N/A              |
| | NRS_AUTOFLAT      |          |                   |                       |                  |
| | NRS_LAMP          |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRS_IMAGE         | N/A      | calwebb_detector1 | calwebb_image2        | N/A              |
| | NRS_WATA          |          |                   |                       |                  |
| | NRS_MSATA         |          |                   |                       |                  |
| | NRS_TACONFIRM     |          |                   |                       |                  |
| | NRS_CONFIRM       |          |                   |                       |                  |
| | NRS_FOCUS         |          |                   |                       |                  |
| | NRS_MIMF          |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRS_FIXEDSLIT     | N/A      | calwebb_detector1 | calwebb_spec2         | calwebb_spec3    |
| | NRS_IFU           |          |                   |                       |                  |
| | NRS_MSASPEC       |          |                   |                       |                  |
+---------------------+----------+-------------------+-----------------------+------------------+
| | NRS_BRIGHTOBJ     | True     | calwebb_tso1      | calwebb_tso-spec2     | calwebb_tso3     |
+---------------------+----------+-------------------+-----------------------+------------------+

Wavefront Sensing and Control Images
------------------------------------
Exposures obtained by any instrument for the purpose of WaveFront Sensing and
Control (WFS&C) use a dedicated processing flow through the pipeline stages.

 - Stage 1: WFS&C exposures use the same :ref:`calwebb_detector1 <calwebb_detector1>`
   pipeline processing and steps as regular images.

 - Stage 2: The ASN generator identifies WFS&C exposures as those with the
   string "WFSC" within their "VISITYPE" keyword value. The resulting ASN
   uses the :ref:`calwebb_wfs-image2 <calwebb_wfs-image2>` pipeline for stage 2
   processing.
   This pipeline configuration is identical to :ref:`calwebb_image2 <calwebb_image2>`
   with the omission of the :ref:`resample <resample_step>` step.

 - Stage 3: The ASN generator identifies pairs of dithered WFS&C images to be
   combined via the "PATTTYPE" keyword value "WFSC". The resulting ASN
   uses the :ref:`calwebb_wfs-image3 <calwebb_wfs-image3>` pipeline
   for stage 3 processing. This pipeline consists of the single step
   :ref:`wfs_combine <wfs_combine_step>`.
