.. _cubepar_reffile:

CUBEPAR reference file
----------------------

:REFTYPE: CUBEPAR
:Data models: `~jwst.datamodels.MiriIFUCubeParsModel`, `~jwst.datamodels.NirspecIFUCubeParsModel`

The CUBEPAR reference file contains parameter values used to construct
the output IFU cubes.

.. include:: ../references_general/cubepar_selection.rst

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
===================  ========  ===========  =============

The formats of the individual table extensions are listed below,
first for the MIRI form of the reference file and then for NIRSpec.

.. include:: ../references_general/cubepar_format.rst

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

