.. _resol_reffile:

RESOL reference file
--------------------

:REFTYPE: RESOL
:Data model: `~jwst.datamodels.MiriResolutionModel`

The RESOL reference file contains the MIRI MRS PSF and LSF widths, per
wavelength band.
This information is used if the ``cube_build`` weight function incorporates the
size of the PSF and LSF, i.e. when using the parameter setting
"--weighting = miripsf".

.. include:: ../references_general/resol_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for RESOL
++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in RESOL reference files,
because they are used as CRDS selectors
(see :ref:`resol_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
CHANNEL    model.meta.instrument.channel
=========  ==============================

Reference File Format
+++++++++++++++++++++
MIRI RESOL reference files are FITS format, with three BINTABLE
extensions.
The FITS primary data array is assumed to be empty.
The format and content of the file is as follows:

===================  ========  ============
EXTNAME              XTENSION  Dimensions
===================  ========  ============
RESOLVING_POWER      BINTABLE  TFIELDS = 11
PSF_FWHM_ALPHA       BINTABLE  TFIELDS = 5
PSF_FWHM_BETA        BINTABLE  TFIELDS = 5
===================  ========  ============

The formats of the individual table extensions are listed below.

.. include:: ../references_general/resol_format.rst

The RESOLVING_POWER table contains information to use for each band.
This table has 12 rows and 11 columns; one row of information for each
MIRI band. The 11 columns contain the polynomial coefficients used to determine
the resolving power for each band. 
The PSF_FWHM_ALPHA table has a format of 1 row and 5 columns. The 5 columns
contain the polynomial coefficients used for determining the alpha coordinate
PSF size. 
Similarly, the PSF_FWHM_BETA table has a format of 1 row and 5 columns,
which contain the polynomial coefficients used for determining the beta 
coordinate PSF size. 
