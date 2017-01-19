Pipeline Classes and Configuration Files
========================================

The actual pipelines that call individual correction steps in various
orders are defined as python classes within python code modules. The pipelines
can be executed by referencing their class name or through the use of a
configuration (.cfg) file that in turn references the class. The table below
shows the pipeline classes that are currently available, the
corresponding pre-defined configurations that make use of those classes, and
the instrument modes to which they can be applied.

+----------------------+------------------------+------------------------------------------+
| Class Name           | Configuration File     | Used For                                 |
+======================+========================+==========================================+
| SloperPipeline       | calwebb_sloper.cfg     | Level-2a processing: all modes           |
+----------------------+------------------------+------------------------------------------+
| DarkPipeline         | calwebb_dark.cfg       | Level-2a processing: darks               |
+----------------------+------------------------+------------------------------------------+
| Image2Pipeline       | calwebb_image2.cfg     | Level-2b processing: imaging modes       |
+----------------------+------------------------+------------------------------------------+
| Spec2Pipeline        | calwebb_spec2.cfg      | Level-2b processing: spectroscopy modes  |
+----------------------+------------------------+------------------------------------------+
| Image3Pipeline       | calwebb_image3.cfg     | Level-3 processing: imaging modes        |
+----------------------+------------------------+------------------------------------------+
| Ami3Pipeline         | calwebb_ami3.cfg       | Level-3 processing: NIRISS AMI mode      |
+----------------------+------------------------+------------------------------------------+

Level-2a Pipelines Step Flow
============================

The list of correction steps applied by the Build 7 calwebb_sloper level-2a (ramps-to-slopes)
pipeline is as follows.

==============  ==============
calwebb_sloper  calwebb_sloper
(All Near-IR)   (MIRI)
==============  ==============
dq_init         dq_init
saturation      saturation
ipc             ipc       
superbias       linearity 
refpix          rscd
linearity       lastframe    
dark_current    dark_current 
\               refpix
jump            jump
ramp_fit        ramp_fit
==============  ==============

Level-2b Imaging Pipeline Step Flow
===================================

The list of correction steps applied by the calwebb_image2 level-2b imaging pipeline
is as follows.

+----------------+
| calwebb_image2 |
+================+
| assign_wcs     |
+----------------+
| flat_field     |
+----------------+
| photom         |
+----------------+


Level-2b Spectroscopic Pipelines Step Flow
==========================================

The list of correction steps invoked by the level-2b spectroscopic
pipeline configurations is shown below. Some steps are only applied to
certain instruments or instrument modes, as noted in the table below.

The calwebb_spec2 pipeline is able to apply background subtraction operations
via either nodded source exposures or dedicated background exposures.
The background-subtraction operations
require the use of associations of exposures. The calwebb_spec2 pipeline
can take either an association table as input, which lists the various 
level-2a products to be processed and their relationships
to one another, or it can take a single level-2a product as input. If an
association table is provided as input, the background-subtraction and
regular calibration steps are applied, in turn, to all input exposures.
If a single level-2a product is used as input, the background subtraction
operations will of course be skipped and only the regular calibration steps
will be applied to the exposure.

The imprint_subtract step, which is only applied to NIRSpec exposures, also
requires the use of an input association table, in order to specify the
imprint exposure to be used in the step. This step will be skipped when
processing a single input product.

The resamp_spec step produces a resampled product for non-IFU modes of
some kinds of spectroscopic exposures. If the resample_spec step is not applied
to a given exposure type, the extract_1d operation will be performed on the
unresampled data.
Cube_build produces a resampled cube for IFU exposures.

Note that level-2b processing for NIRCam and NIRISS Wide-Field Slitless (grism)
Spectroscopy modes is not yet implemented.

+------------------+----+-----+-----+----+----+-----+--------+
| Instrument Mode  |     NIRSpec    |     MIRI      | NIRISS |
+------------------+----+-----+-----+----+----+-----+--------+
| Step             | FS | MOS | IFU | FS | SL | MRS |  SOSS  |
+==================+====+=====+=====+====+====+=====+========+
| assign_wcs       | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| bkg_subtract     | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| imprint_subtract |    |  X  |  X  |    |    |     |        |
+------------------+----+-----+-----+----+----+-----+--------+
| extract_2d       | X  |  X  |     |    |    |     |        |
+------------------+----+-----+-----+----+----+-----+--------+
| flat_field       | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| srctype          | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| straylight       |    |     |     |    |    |  X  |        |
+------------------+----+-----+-----+----+----+-----+--------+
| fringe           |    |     |     |    |    |  X  |        |
+------------------+----+-----+-----+----+----+-----+--------+
| photom           | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+
| resamp_spec      | X  |  X  |     | X  |    |     |        |
+------------------+----+-----+-----+----+----+-----+--------+
| cube_build       |    |     |  X  |    |    |  X  |        |
+------------------+----+-----+-----+----+----+-----+--------+
| extract_1d       | X  |  X  |  X  | X  | X  |  X  |   X    |
+------------------+----+-----+-----+----+----+-----+--------+

Level-3 Imaging Pipeline Step Flow
==================================

The correction steps applied by the calwebb_image3 pipeline are shown
below. The tweakreg_catalog, tweakreg, skymatch, and outlier_detection steps
can only be applied when an association table containing multiple inputs
is given. If a single input product is specified, only the resample and
source_catalog steps will be applied.

+-------------------+
| calwebb_image3    |
+===================+
| tweakreg_catalog  |
+-------------------+
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

Level-3 Aperture Masking Interferometry (AMI) Pipeline Step Flow
================================================================

The correction steps applied by the calwebb_ami3 pipeline are shown
below.

+---------------+
| calwebb_ami3  |
+===============+
| ami_analyze   |
+---------------+
| ami_average   |
+---------------+
| ami_normalize |
+---------------+

