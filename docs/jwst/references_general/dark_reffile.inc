.. _dark_reffile:

DARK Reference File
-------------------

:REFTYPE: DARK
:Data models: `~jwst.datamodels.DarkModel`, `~jwst.datamodels.DarkMIRIModel`

The DARK reference file contains pixel-by-pixel and frame-by-frame
dark current values for a given detector readout mode.

.. include:: ../references_general/dark_selection.inc

.. include:: ../includes/standard_keywords.inc

.. include:: ../references_general/dark_specific.inc
             
Reference File Format
+++++++++++++++++++++
DARK reference files are FITS format, with 3 IMAGE extensions
and 1 BINTABLE extension. The FITS primary data array is assumed to be empty.
The format and content of the files is different for MIRI than the
near-IR instruments, as shown below.

Near-IR Detectors
~~~~~~~~~~~~~~~~~
Characteristics of the three IMAGE extensions for DARK files used with
the Near-IR instruments are as follows (see `~jwst.datamodels.DarkModel`):

=======  =====  =======================  =========
EXTNAME  NAXIS  Dimensions               Data type
=======  =====  =======================  =========
SCI      3      ncols x nrows x ngroups  float
ERR      3      ncols x nrows x ngroups  float
DQ       2      ncols x nrows            integer
DQ_DEF   2      TFIELDS = 4              N/A
=======  =====  =======================  =========

MIRI Detectors
~~~~~~~~~~~~~~
The DARK reference files for the MIRI detectors depend on the integration number,  
because the first integration of MIRI exposures contains effects from the detector
reset and are slightly different from subsequent integrations. Currently the MIRI
DARK reference files contain a correction for only two integrations: the first
integration of the DARK is subtracted from the first integration of the science data,
while the second DARK integration is subtracted from all subsequent science integrations.
The format of the MIRI DARK reference files is as follows
(see `~jwst.datamodels.DarkMIRIModel`):

=======  =====  ===============================  =========
EXTNAME  NAXIS  Dimensions                       Data type
=======  =====  ===============================  =========
SCI      4      ncols x nrows x ngroups x nints  float
ERR      4      ncols x nrows x ngroups x nints  float
DQ       4      ncols x nrows x 1 x nints        integer
DQ_DEF   2      TFIELDS = 4                      N/A
=======  =====  ===============================  =========

.. include:: ../includes/dq_def.inc
