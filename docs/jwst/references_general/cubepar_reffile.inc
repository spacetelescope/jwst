.. _cubepar_reffile:

CUBEPAR reference file
----------------------

:REFTYPE: CUBEPAR
:Data models: `~jwst.datamodels.MiriIFUCubeParsModel`, `~jwst.datamodels.NirspecIFUCubeParsModel`

The CUBEPAR reference file contains parameter values used to construct
the output IFU cubes.

.. include:: ../references_general/cubepar_selection.inc

.. include:: ../includes/standard_keywords.inc

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

MIRI Reference File Format
++++++++++++++++++++++++++
The MIRI CUBEPAR reference files are FITS format, with  5  BINTABLE extensions.
The FITS primary data array is assumed to be empty.
The format and content of the MIRI CUBEPAR reference file 

===================  ========  ===========  
EXTNAME              XTENSION  Dimensions   
===================  ========  ===========  
CUBEPAR              BINTABLE  TFIELDS = 6  
CUBEPAR_MSM          BINTABLE  TFIELDS = 6  
MULTICHANNEL_MSM     BINTABLE  TFIELDS = 5  
CUBEPAR_EMSM         BINTABLE  TFIELDS = 5  
MULTICHANNEL_EMSM    BINTABLE  TFIELDS = 4
MULTICHANNEL_DRIZ    BINTABLE  TFIELDS = 1
===================  ========  ===========  


NIRSPec Reference File Format
+++++++++++++++++++++++++++++
The NIRSpec CUBEPAR reference files are FITS format, with 9  BINTABLE extensions.

====================  ========  ===========
EXTNAME               XTENSION  Dimensions 
====================  ========  ===========
CUBEPAR               BINTABLE  TFIELDS = 6
CUBEPAR_MSM           BINTABLE  TFIELDS = 6
MULTICHAN_PRISM_MSM   BINTABLE  TFIELDS = 5
MULTICHAN_MED_MSM     BINTABLE  TFIELDS = 5
MULTICHAN_HIGH_MSM    BINTABLE  TFIELDS = 5
CUBEPAR_EMSM          BINTABLE  TFIELDS = 5
MULTICHAN_PRISM_EMSM  BINTABLE  TFIELDS = 4
MULTICHAN_MED_EMSM    BINTABLE  TFIELDS = 4
MULTICHAN_HIGH_EMSM   BINTABLE  TFIELDS = 4
====================  ========  ===========

The formats of the individual table extensions are listed below,
first for the MIRI  reference file and then for NIRSpec.

.. include:: ../references_general/cubepar_format.inc

These reference files contain tables for each wavelength band giving the spatial
and spectral size, and the size of the region of interest (ROI) to use to
construct an IFU cube.  If only one band is used to construct the IFU cube then the *CUBEPAR* and *CUBEPAR_MSM* or *CUBE_EMSM* tables are used. These types of cubes will have a linear - wavelength dimension.  
If more than one wavelength band is used to build the IFU cube then the MULTICHANNEL (MIRI)
or MULTICHAN (NIRSPEC) tables are used o set the spectral and spatial roi size, and the wavelength dependent weighting function parameters.
For multi-band IFU cubes  then the final
spatial size will be the smallest one from the list of input bands and these cubes will have a non-linear wavelength dimension. 


The MIRI reference table descriptions:

- **CUBEPAR** table contains the spatial and spectral cube sample size for each band.
- **CUBEPAR_MSM** table contains the Modified Shepard Method (MSM) weighting values to use for each band.
- **MULTICHANNEL_MSM** table is used for the MSM weighting and  contains the wavelengths and associated region of interest size to use when IFU cubes are created from several bands and the final output is to have an IFU cube of varying spectral scale.

- **CUBEPAR_EMSM** table contains the Exponential Modified Shepard Method (EMSM) weighting values to use for each band.
- **MULTICHANNEL_EMSM** table is used for the EMSM weighting and contains the wavelengths and associated region of interest size to use when IFU cubes are created from several bands and the final output is to have an IFU cube of varying spectral scale.
- **MULTICHANNEL_DRIZ** table is used for the DRIZZLE weighting and contains the wavelengths to use when IFU cubes are created from several bands and the final output is to have an IFU cube of varying spectral scale.


The NIRSPEC reference table descriptions:

- **CUBEPAR** table contains the spatial and spectral cube sample size for each band.
- **CUBEPAR_MSM** table contains the Modified Shepard Method (MSM) weighting values to use for each band.

- **MULTICHAN_PRISM_MSM** table is used for the MSM weighting and  contains the wavelengths and associated region of interest size to use when IFU cubes are created from the grating prism and the final IFU Cube output has a varying spectral scale.

- **MULTICHAN_MED_MSM** table is used for the MSM weighting and  contains the wavelengths and associated region of interest size to use when IFU cubes are created from the medium resolution grating and the final IFU Cube output has a varying spectral scale.

- **MULTICHAN_HIGH_MSM** table is used for the MSM weighting and  contains the wavelengths and associated region of interest size to use when IFU cubes are created from the high resolution gratings and the final IFU Cube output has a varying spectral scale.

- **CUBEPAR_EMSM** table contains the Exponential Modified Shepard Method (EMSM) weighting values to use for each band.

- **MULTICHAN_PRISM_EMSM** table is used for the EMSM weighting and  contains the wavelengths and associated region of interest size to use when IFU cubes are created from the grating prism and the final IFU Cube output has a varying spectral scale.

- **MULTICHAN_MED_EMSM** table is used for the EMSM weighting and  contains the wavelengths and associated region of interest size to use when IFU cubes are created from the medium resolution grating and the final IFU Cube output has a varying spectral scale.

- **MULTICHAN_HIGH_EMSM** table is used for the EMSM weighting and  contains the wavelengths and associated region of interest size to use when IFU cubes are created from the high resolution gratings and the final IFU Cube output has a varying spectral scale.

