.. _mrsptcorr_reffile:

MRSPTCORR reference file
--------------------------

:REFTYPE: MRSPTCORR
:Data models: `~jwst.datamodels.MirMrsPtCorrModel`

The MRSPTCORR reference file contains parameter values used to subtract
the MRS 12 micron spectral leak in the spectral leak step. It also contains parameters
to correct point sources, in future enhancements,  for  across-slice corrections and
throughput variations. 

.. include:: ../references_general/mrsptcorr_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for MRSPTCORR
++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in MRSPTCORR reference files,
because they are used as CRDS selectors
(see :ref:`mrsptcorr_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

MIRI Reference File Format
++++++++++++++++++++++++++
The MIRI MRSPTCORR reference files are FITS format, with  5  BINTABLE extensions.
The FITS primary data array is assumed to be empty.
The format and content of the MIRI MRSPTCORR reference file

===================  ========  ===========  
EXTNAME              XTENSION  Dimensions   
===================  ========  ===========  
LEAKCOR              BINTABLE  TFIELDS = 3
TRACOR               BINTABLE  TFIELDS = 7
WAVCORR_OPTICAL      BINTABLE  TFIELDS = 6
WAVCORR_XSLICE       BINTABLE  TFIELDS = 2
WAVCORR_SHIFT        BINTABLE  TFIELDS = 3
===================  ========  ===========  

The formats of the individual table extensions are listed below.

.. include:: ../references_general/mrsptcorr_format.inc

This reference file contains the relative spectrophotometric response function for the MRS 12 micron
spectral leak, along with tables of across-slice wavelength and throughput variations for point
sources in each of the MRS bands. Currently, only the LEAKCOR table is use in the jwst pipeline. 

