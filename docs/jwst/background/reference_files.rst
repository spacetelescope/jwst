Reference Files
===============
The background image subtraction step uses reference files only when
processing Wide-Field Slitless Spectroscopy (WFSS) exposures. Two reference
files are used for WFSS mode.

WFSS Background reference file
------------------------------

:REFTYPE: WFSSBKG
:Data model: `WfssBkgModel`

The WFSS background reference file contains a "master" image of the
dispersed background produced by a particular filter+grism combination.

.. include:: selection_criteria.rst
 
.. include:: ../includes/standard_keywords.rst


Type Specific Keywords for WFSSBKG
++++++++++++++++++++++++++++++++++
The following additional keywords are required for the WFSSBKG reference
type:

=========  ========================
Keyword    Model Name
=========  ========================
DETECTOR   meta.instrument.detector
EXP_TYPE   meta.exposure.type
FILTER     meta.instrument.filter
PUPIL      meta.instrument.pupil
=========  ========================

Reference File Format
+++++++++++++++++++++
WFSSBKG reference files are FITS files with 3 IMAGE extensions and
1 BINTABLE extension. The FITS primary data array is assumed to be empty.
The characteristics of the FITS extensions are as follows:

=======  ========  =====  ==============  =========
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  =========
SCI      IMAGE       2    ncols x nrows   float
ERR      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   integer
DQ_DEF   BINTABLE    2    TFIELDS = 4     N/A
=======  ========  =====  ==============  =========

.. include:: ../includes/dq_def.rst

Wavelength Range reference file
-------------------------------

:REFTYPE: WAVELENGTHRANGE
:Data model: `WavelengthrangeModel`

The wavelength range reference file contains information about the range of
wavelengths in the exposure. It is used, together with a source catalog,
to create a mask giving the locations of source spectra in the target image
and hence where the background regions are.

CRDS Selection Criteria
+++++++++++++++++++++++
Wavelengthrange reference files are selected by:

  INSTRUME, EXP_TYPE, PUPIL (NIRCam only), and MODULE (NIRCam only)

