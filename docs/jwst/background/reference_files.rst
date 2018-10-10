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

.. include:: wfssbkg_selection.rst
 
.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for WFSSBKG
++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in WFSSBKG reference files,
because they are used as CRDS selectors
(see :ref:`wfssbkg_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
EXP_TYPE   model.meta.exposure.type
FILTER     model.meta.instrument.filter
PUPIL      model.meta.instrument.pupil
=========  ==============================

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

.. include:: wavelengthrange_selection.rst

Standard Keywords
+++++++++++++++++
**NOTE:** WAVELENGTHRANGE requires the standard keywords shown above.


Type Specific Keywords for WAVELENGTHRANGE
++++++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following additional keywords are *required* in WAVELENGTHRANGE
reference files, because they are used as CRDS selectors
(see :ref:`wavelengthrange_selectors`):

=========  ========================
Keyword    Data Model Name
=========  ========================
EXP_TYPE   model.meta.exposure.type
=========  ========================

