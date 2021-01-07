.. _area_reffile:

AREA Reference File
-------------------

:REFTYPE: AREA

The AREA reference file contains pixel area information for a given
instrument mode.

.. include:: ../references_general/area_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for AREA
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in AREA reference files,
because they are used as CRDS selectors
(see :ref:`area_selectors`):

=========  ==============================  =============================
Keyword    Data Model Name                 Instrument
=========  ==============================  =============================
DETECTOR   model.meta.instrument.detector  All
EXP_TYPE   model.meta.exposure.type        MIRI, NIRCam, NIRISS, NIRSpec
FILTER     model.meta.instrument.filter    MIRI, NIRCam, NIRISS, NIRSpec
PUPIL      model.meta.instrument.pupil     NIRCam, NIRISS
GRATING    model.meta.instrument.grating   NIRSpec
=========  ==============================  =============================

Reference File Format
+++++++++++++++++++++
AREA reference files are FITS format.
For imaging modes (FGS, MIRI, NIRCam, and NIRISS) the AREA reference files
contain 1 IMAGE extension, while reference files for NIRSpec spectroscopic
modes contain 1 BINTABLE extension.
The FITS primary HDU does not contain a data array.

Imaging Modes
^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.PixelAreaModel`

The format of imaging mode AREA reference files is as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
=======  ========  =====  ==============  =========

The SCI extension data array contains a 2-D pixel-by-pixel map of relative
pixel areas, normalized to a value of 1.0. The absolute value of the nominal
pixel area is given in the primary header keywords PIXAR_SR and PIXAR_A2, in
units of steradians and square arcseconds, respectively.
These keywords should have the same values as the
PIXAR_SR and PIXAR_A2 keywords in the header of the corresponding PHOTOM
reference file.

NIRSpec Fixed-Slit Mode
^^^^^^^^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.NirspecSlitAreaModel`

The BINTABLE extension has EXTNAME='AREA' and has column characteristics
shown below. There is one row for each of the 5 fixed slits, with ``slit_id``
values of "S200A1", "S200A2", "S400A1", "S200B1", and "S1600A1". The pixel
area values are in units of square arcseconds and represent the nominal
area of any pixel illuminated by the slit.

===========  =========
Column name  Data type 
===========  =========
slit_id      char\*15
pixarea      float
===========  =========

NIRSpec MOS Mode
^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.NirspecMosAreaModel`

The BINTABLE extension has EXTNAME='AREA' and has column characteristics
shown below. There is one row for each shutter in each MSA quadrant. The
quadrant and shutter values are 1-indexed. The pixel area values are in
units of square arcseconds and represent the nominal area of any pixel
illuminated by a given MSA shutter.

===========  =========
Column name  Data type 
===========  =========
quadrant     integer
shutter_x    integer
shutter_y    integer
pixarea      float
===========  =========

NIRSpec IFU Mode
^^^^^^^^^^^^^^^^

:Data model: `~jwst.datamodels.NirspecIfuAreaModel`

The BINTABLE extension has EXTNAME='AREA' and has column characteristics
shown below. There is one row for each of the 30 IFU slices, with the
``slice_id`` values being 0-indexed (i.e. range from 0 to 29). The pixel
area values are in units of square arcseconds and represent the nominal
area of any pixel illuminated by a given IFU slice.

===========  =========
Column name  Data type 
===========  =========
slice_id     integer
pixarea      float
===========  =========

