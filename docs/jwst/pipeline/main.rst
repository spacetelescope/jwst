.. _pipelines:

Pipeline Modules
================

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

In addition, there are several pipeline modules designed for special processing,
including:

- :ref:`calwebb_dark` for processing dark exposures
- :ref:`calwebb_guider` for processing FGS guide star data

Each pipeline consists of a certain sequence of calibration steps and is
defined as a python class within a python code module. The pipelines
can be executed from the command line either by referencing their class name or
by supplying a configuration (.cfg) file that in turn references the pipeline class.
From within python, the pipelines are called by their class names, but
configuration files can still be supplied in order to set pipeline or step
parameter values.
The table below lists the pipeline classes that are currently available, the
corresponding configuration files that call those classes, and
the observing modes for which they are intended.

+-------------------+---------------------------+------------------------------+
| Class Name        | Configuration File        | Used For                     |
+===================+===========================+==============================+
| Detector1Pipeline | calwebb_detector1.cfg     | Stage 1: all non-TSO modes   |
+-------------------+---------------------------+------------------------------+
| Detector1Pipeline | calwebb_tso1.cfg          | Stage 1: all TSO modes       |
+-------------------+---------------------------+------------------------------+
| DarkPipeline      | calwebb_dark.cfg          | Stage 1: darks               |
+-------------------+---------------------------+------------------------------+
| GuiderPipeline    | calwebb_guider.cfg        | Stage 1+2: FGS guiding modes |
+-------------------+---------------------------+------------------------------+
| Image2Pipeline    | calwebb_image2.cfg        | Stage 2: imaging modes       |
+-------------------+---------------------------+------------------------------+
| Image2Pipeline    | calwebb_tso-image2.cfg    | Stage 2: TSO imaging modes   |
+-------------------+---------------------------+------------------------------+
| Image2Pipeline    | calwebb_wfs-image2.cfg    | Stage 2: WFS&C imaging       |
+-------------------+---------------------------+------------------------------+
| Spec2Pipeline     | calwebb_spec2.cfg         | Stage 2: spectroscopy modes  |
+-------------------+---------------------------+------------------------------+
| Spec2Pipeline     | calwebb_tso-spec2.cfg     | Stage 2: TSO spectral modes  |
+-------------------+---------------------------+------------------------------+
| Spec2Pipeline     | calwebb_nrslamp-spec2.cfg | Stage 2: NIRSpec lamps       |
+-------------------+---------------------------+------------------------------+
| Image3Pipeline    | calwebb_image3.cfg        | Stage 3: imaging modes       |
+-------------------+---------------------------+------------------------------+
| WfsCombineStep    | calwebb_wfs-image3.cfg    | Stage 3: WFS&C imaging       |
+-------------------+---------------------------+------------------------------+
| Spec3Pipeline     | calwebb_spec3.cfg         | Stage 3: spectroscopy modes  |
+-------------------+---------------------------+------------------------------+
| Ami3Pipeline      | calwebb_ami3.cfg          | Stage 3: NIRISS AMI mode     |
+-------------------+---------------------------+------------------------------+
| Coron3Pipeline    | calwebb_coron3.cfg        | Stage 3: Coronagraphic mode  |
+-------------------+---------------------------+------------------------------+
| TSO3Pipeline      | calwebb_tso3.cfg          | Stage 3: TSO modes           |
+-------------------+---------------------------+------------------------------+

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

