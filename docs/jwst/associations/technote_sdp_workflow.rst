.. _sdp-workflow:

================================
Science Data Processing Workflow
================================

General Workflow
================

See :ref:`level3-asn-jwst-overview` for an overview of how JWST uses
associations. This document describes how associations are used by the
ground processing system to execute the Level 2 and Level 3 pipelines
based on.

Up to the initial calibration step `calwebb_sloper`, the science
exposures are treated individually. However, starting at the Level 2
calibration step, exposures may need other exposures in order to be
further processed. Instead of creating a single monolithic pipeline,
the workflow uses the associations to determine what pipeline should
be executed and when to execute them. In the figure below, this
wait-then-execute process is represented by the `workflow trigger`.
The workflow reads the contents of an association to determine what
exposures, and possibly other files, are needed to continue
processing. The workflow then waits until all exposures exist. At that
point, the related calibration step is executed with the association
as input.

With this finer granularity, the workflow can run more processes parallel,
and allows the operators deeper visibility into the progression of the
data through the system.

.. figure:: graphics/workflow-generic.png
   :scale: 75%

   General workflow through Level 2 and Level 3 processing

The figure represents the following workflow:

- Data comes down from the observatory and is converted to the Level
  1a FITS files.
- `calwebb_sloper` is run on each file to convert the data to the
  Level 2a format.
- In parallel with `calwebb_sloper`, the Pool Maker collects the list
  of downloaded exposures and places them in the Association Pool
- When enough exposures have been download to complete and Association
  Candidate, such as an Observation Candidate, the Pool Maker calls
  the Association Generator, `asn_generate`, to create the set of
  associations based on that Candidate.
- For each association generated, the workflow creates a file watch
  list from the association, then waits until all exposures needed by
  that association come into existence.
- When all exposures for an association exist, the workflow then
  executes the corresponding pipeline, passing the association as
  input.

Wide Field Slitless Spectroscopy
================================

In most cases, the data will flow from Level 2 to Level 3, completing
calibration. However, more complicated situations can be handled by
the same wait-then-execute process. One particular case is for the
Wide Field Slitless Spectrometry (WFSS) modes. The specific flow is
show in the figure below:

.. figure:: graphics/workflow-wfss.png
   :scale: 75%

   WFSS workflow

For WFSS data, at least two observations are made, one consisting of a
direct image of the field-of-view (FOV), and a second where the FOV is
dispersed using a grism. The direct image is first processed through
Level 3. At the Level 3 stage, a source catalog of objects found in
the image, and a segment map, are generated. These files are then used
as input to the Level 2 processing of the spectral data. This extra
link between the two major stages is represented by the `Segment &
Catalog` file set, show in red in the diagram. The level 2 association
`grism_spec2_asn` not only lists the needed level 2a exposures, but
also the catalog and segment map files produced by the Level 3 image
processing. Hence, the workflow knows to wait for these files before
continuing the spectral processing.
