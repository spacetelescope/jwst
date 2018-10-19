Reference Files
===============

The ``cube_build`` step uses two reference files: CUBEPAR and RESOL.
The RESOL reference file is only used for processing MIRI IFU data.
CUBEPAR is used for both NIRSpec and MIRI IFU data.

CUBEPAR reference file
----------------------

:RETYPE: CUBEPAR
:Data models: `~jwst.datamodels.MiriIFUCubeParsModel`, `~jwst.datamodels.NirspecIFUCubeParsModel`

The CUBEPAR reference file contains parameter values used to construct
the output IFU cubes.

.. include:: cubepar_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for CUBEPAR
+++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in CUBEPAR reference files,
because they are used as CRDS selectors
(see :ref:`cubepar_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
+++++++++++++++++++++
CUBEPAR reference files are FITS format, with either 3 (MIRI) or
5 (NIRSpec) BINTABLE extensions.
The FITS primary data array is assumed to be empty.
The format and content of the file is as follows:

===================  ========  ===========  =============
EXTNAME              XTENSION  Dimensions   Instrument
===================  ========  ===========  =============
CUBEPAR              BINTABLE  TFIELDS = 6  Both
CUBEPAR_MSM          BINTABLE  TFIELDS = 6  Both
MULTICHANNEL_MSM     BINTABLE  TFIELDS = 5  MIRI only
MULTICHAN_PRISM_MSM  BINTABLE  TFIELDS = 5  NIRSpec only
MULTICHAN_MED_MSM    BINTABLE  TFIELDS = 5  NIRSpec only
MULTICHAN_HIGH_MSM   BINTABLE  TFIELDS = 5  NIRSpec only
===================  ========  ===========  ============

The formats of the individual table extensions are listed below,
first for the MIRI form of the reference file and then for NIRSpec.

.. include:: cubepar_format.rst

These reference files contain tables for each wavelength band giving the spatial
and spectral size, and the size of the region of interest (ROI) to use to
construct an IFU cube.
If more than one wavelength band is used to build the IFU cube, then the final
spatial and spectral size will be the smallest one from the list of input bands.
The "CUBEPAR" table contains the spatial and spectral cube sample size for each
band. The "CUBEPAR_MSM" table contains the Modified Shepard Method (MSM)
weighting values to use for each band.
The "MULTICHANNEL_MSM" table contains the wavelengths and associated
region of interest size to use when IFU cubes are created from several bands
and the final output is to have an IFU cube of varying spectral scale.

For MIRI, the twelve spectral bands can be combined into a single IFU cube,
in which case all of the information needed to create cubes of varying
wavelength sampling is contained in this table.
For NIRSpec IFU data, however, there are three types of multi-band cubes that
can be created: PRISM, MEDIUM, and HIGH resolution.
The three MULTICHAN_<type>_MSM tables in the NIRSpec version of the reference
file contain the wavelength sampling and region of interest size
information to use with each of these types of multi-band cubes.

RESOL reference file
--------------------

:RETYPE: RESOL
:Data model: `~jwst.datamodels.MiriResolutionModel`

The RESOL reference file contains the MIRI MRS PSF and LSF widths, per
wavelength band.
This information is used if the ``cube_build`` weight function incorporates the
size of the PSF and LSF, i.e. when using the parameter setting
"--weighting = miripsf".

.. include:: resol_selection.rst

.. include:: ../includes/standard_keywords.rst

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

.. include:: resol_format.rst

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
