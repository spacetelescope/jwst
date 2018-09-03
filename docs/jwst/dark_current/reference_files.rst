Reference File
==============
The dark current step uses a DARK reference file.

.. include:: dark_selection.rst
             
.. include:: ../includes/standard_keywords.rst

.. include:: dark_specific.rst
             
DARK Reference File Format
--------------------------

NIR Detectors
+++++++++++++

Dark reference files are FITS files with 3 IMAGE extensions and 1 BINTABLE
extension. The FITS primary data array is assumed to be empty. The 
characteristics of the three image extensions for darks used with the
Near-IR instruments are as follows:

=======  =====  =======================  =========
EXTNAME  NAXIS  Dimensions               Data type
=======  =====  =======================  =========
SCI      3      ncols x nrows x ngroups  float
ERR      3      ncols x nrows x ngroups  float
DQ       2      ncols x nrows            integer
DQ_DEF   2      TFIELDS = 4              N/A
=======  =====  =======================  =========

.. include:: ../includes/dq_def.rst

MIRI Detectors
++++++++++++++

The dark reference files for the MIRI detectors depend on the integration number,  
because the first integration of MIRI exposures contains effects from the detector
reset and are slightly different from subsequent integrations. Currently the MIRI
dark reference files contain a correction for only two integrations: the first
integration of the dark is subtracted from the first integration of the science data,
while the second dark integration is subtracted from all subsequent science integrations.
The format of the MIRI dark reference files is as follows:

=======  =====  ===============================  =========
EXTNAME  NAXIS  Dimensions                       Data type
=======  =====  ===============================  =========
SCI      4      ncols x nrows x ngroups x nints  float
ERR      4      ncols x nrows x ngroups x nints  float
DQ       4      ncols x nrows x 1 x nints        integer
DQ_DEF   2      TFIELDS = 4                      N/A
=======  =====  ===============================  =========

.. include:: ../includes/dq_def.rst
