.. _fringefreq_reffile:

FRINGEFREQ reference file
-------------------------

:REFTYPE: FRINGEFREQ
:Data models: `~jwst.datamodels.FringeFreqModel`

The FRINGEFREQ reference files contain parameter values used to correct MIRI MRS images
for  residual fringes that remain  after applying the fringe flat.

.. include:: ../references_general/fringefreq_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for FRINGEFREQ
+++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in FRINGEFREQ reference files,
because they are used as CRDS selectors
(see :ref:`fringefreq_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
++++++++++++++++++++++++++
The FRINGEFREQ reference files are FITS format, with  4  BINTABLE extensions.
The FITS primary data array is assumed to be empty.
The format and content of the FRINGEFREQ reference file is

===================  ========  ===========  
EXTNAME              XTENSION  Dimensions   
===================  ========  ===========  
RFC_FREQ_SHORT       BINTABLE  TFIELDS = 7  
RFC_FREQ_MEDIUM      BINTABLE  TFIELDS = 7  
RFC_FREQ_LONG        BINTABLE  TFIELDS = 7  
MAX_AMP              BINTABLE  TFIELDS = 2  
===================  ========  ===========  

The formats of the individual table extensions are listed below.

.. include:: ../references_general/fringefreq_format.inc

These reference files contain tables for each MIRI band giving the
fringe frequencies and other parameters for each band to fit and remove residual fringes. 


The reference table descriptions:

- **RFC_FREQ_SHORT** table contains the fringe frequencies and parameters for the SHORT band.
- **RFC_FREQ_MEDIUM** table contains the fringe frequencies and parameters for the MEDIUM band.
- **RFC_FREQ_LONG** table contains the fringe frequencies and parameters for the LONG band.
- **MAX_AMP** table contains a wavelength dependent maximum amplitude which is use for feature identification and fit rejection. 



