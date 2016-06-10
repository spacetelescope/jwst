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
| Spec2Pipeline        | calniriss_soss2.cfg    | Level-2b processing: NIRISS SOSS         |
+----------------------+------------------------+------------------------------------------+
| Image3Pipeline       | calwebb_image3.cfg     | Level-3 processing: imaging modes        |
+----------------------+------------------------+------------------------------------------+

Level-2a Pipelines Step Flow
============================

The list of correction steps invoked by the Build 6 level-2a (ramps-to-slopes)
pipeline configurations is as follows.

==============  ==============
calwebb_sloper  calwebb_sloper
(All Near-IR)   (MIRI)
==============  ==============
dq_init         dq_init
saturation      saturation
ipc             ipc       
superbias       \         
refpix          refpix
\               reset    
\               lastframe    
linearity       linearity    
dark_current    dark_current 
jump            jump
ramp_fit        ramp_fit
==============  ==============

Level-2b Imaging Pipeline Step Flow
===================================

The list of correction steps invoked by the level-2b imaging pipeline
configuration is as follows. Note that the persistence and emission
steps are currently no-op (no operation applied to the data) placeholders.

+----------------+
| calwebb_image2 |
+================+
| assign_wcs     |
+----------------+
| flat_field     |
+----------------+
| persistence    |
+----------------+
| emission       |
+----------------+
| photom         |
+----------------+


Level-2b Spectroscopic Pipelines Step Flow
==========================================

The list of correction steps invoked by the level-2b spectroscopic
pipeline configurations is shown below. The calniriss_soss2 configuration
skips the assign_wcs step, because the WCS for SOSS mode is not
implemented. This means that the extract_2d step is also currently
skipped for NIRISS SOSS exposures, because that operation requires WCS
information.

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

The straylight and fringe steps are only applied to MIRI MRS exposures.

================ ===============
calwebb_spec2    calniriss_soss2
================ ===============
assign_wcs       \
bkg_subtract     bkg_subtract
imprint_subtract \
extract_2d       \
flat_field       flat_field
straylight       \
fringe           \
photom           photom
================ ===============

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

