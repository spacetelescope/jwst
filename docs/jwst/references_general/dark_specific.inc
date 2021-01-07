Type Specific Keywords for DARK
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in DARK reference files,
because they are used as CRDS selectors
(see :ref:`dark_selectors`):

=========  ==============================  ==========================
Keyword    Data Model Name                 Instruments
=========  ==============================  ==========================
DETECTOR   model.meta.instrument.detector  All
READPATT   model.meta.exposure.readpatt    FGS, MIRI, NIRISS, NIRSpec
SUBARRAY   model.meta.subarray.name        All
SUBSTRT1   model.meta.subarray.xstart      NIRSpec
SUBSTRT2   model.meta.subarray.ystart      NIRSpec
SUBSIZE1   model.meta.subarray.xsize       NIRSpec
SUBSIZE2   model.meta.subarray.ysize       NIRSpec
=========  ==============================  ==========================

