.. _fringefreq_reffile:

FRINGEFREQ reference file
----------------------

:REFTYPE: FRINGEFREQ
:Data models: `~jwst.datamodels.MiriIFUCubeParsModel`

The FRINGEFREQ reference file contains parameter values used to correct MIRI MRS images
for a residual fringe that remains  after applying the fringe flat.

.. include:: ../references_general/fringefreq_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for FRINGEFREQ
+++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in FRINGEFREQ reference files,
because they are used as CRDS selectors
(see :ref:`cubepar_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
++++++++++++++++++++++++++
The FRINGEFREQ reference files are FITS format, with  $  BINTABLE extensions.
The FITS primary data array is assumed to be empty.
The format and content of the MIRI CUBEPAR reference file 

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

These reference files contain tables for each MIRI band  giving the



The reference table descriptions:

- **RFC_FREQ_SHORT** table contains the spatial and spectral cube sample size for each band.
- **RFC_FREQ_MEDIUM** table contains the Modified Shepard Method (MSM) weighting values to use for each band.
- **RFC_FREQ_LONG** table is used for the MSM weighting and  contains the wavelengths and associated region of interest size to use when IFU cubes are created from several bands and the final output is to have an IFU cube of varying spectral scale.

- **MAX_AMP** table contains the Exponential Modified Shepard Method (EMSM) weighting values to use for each band.



