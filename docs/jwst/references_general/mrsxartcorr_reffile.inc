.. _mrsxartcorr_reffile:

MRSXARTCORR reference file
--------------------------

:REFTYPE: MRSXARTCORR
:Data models: `~jwst.datamodels.MirMrsXArtCorrModel`

The MRSXARTCORR reference file contains parameter values used to model and subtract
the cross-artifact in the straylight step.

.. include:: ../references_general/mrsxartcorr_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for MRSXARTCORR
++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in MRSXARTCORR reference files,
because they are used as CRDS selectors
(see :ref:`mrsxartcorr_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

MIRI Reference File Format
++++++++++++++++++++++++++
The MIRI MRSXARTCORR reference files are FITS format, with  12  BINTABLE extensions.
The FITS primary data array is assumed to be empty.
The format and content of the MIRI MRSXARTCORR reference file

===================  ========  ===========  
EXTNAME              XTENSION  Dimensions   
===================  ========  ===========  
1A                   BINTABLE  TFIELDS = 7
1B                   BINTABLE  TFIELDS = 7
1C                   BINTABLE  TFIELDS = 7
2A                   BINTABLE  TFIELDS = 7
2B                   BINTABLE  TFIELDS = 7
2C                   BINTABLE  TFIELDS = 7
3A                   BINTABLE  TFIELDS = 7
3B                   BINTABLE  TFIELDS = 7
3C                   BINTABLE  TFIELDS = 7
4A                   BINTABLE  TFIELDS = 7
4B                   BINTABLE  TFIELDS = 7
4C                   BINTABLE  TFIELDS = 7
===================  ========  ===========  

The formats of the individual table extensions are listed below.

.. include:: ../references_general/mrsxartcorr_format.inc

These reference files contain tables for each wavelength band giving the appropriate kernel properties
to use to model the cross-artifact for each band.  These include Lorentzian plus Gaussian models.

The MIRI reference table contains parameters for each band in the corresponding EXTNAME extension.
