:orphan:

.. _specwcs_reffile:
  
SPECWCS Reference File
----------------------

:REFTYPE: SPECWCS
:Data model: `~jwst.datamodels.SpecwcsModel`

Reference Selection Keywords for SPECWCS
++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate SPECWCS references based on the following keywords.
SPECWCS is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== =========================================================================
Instrument Keywords
========== =========================================================================
MIRI       INSTRUME, EXP_TYPE, DETECTOR, CHANNEL, BAND, SUBARRAY, DATE-OBS, TIME-OBS
NIRCam     INSTRUME, EXP_TYPE, MODULE, PUPIL, DATE-OBS, TIME-OBS
NIRISS     INSTRUME, EXP_TYPE, FILTER, PUPIL, DATE-OBS, TIME-OBS
========== =========================================================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The exact format and contents of the SPECWCS reference file varies depending on the
instrument observing mode, as explained in the following sections.

MIRI LRS mode
:::::::::::::
For the MIRI LRS mode the SPECWCS file is in FITS format.
The reference file contains the zero point offset for the slit relative to the full field of view.
For the Fixed Slit exposure type the zero points in X and Y are stored in the header of the second HDU in the
'IMX' and 'IMY' keywords. For the Slitless exposure type they are stored in the header of the second HDU in
FITS keywords 'IMXSLTl' and 'IMYSLTl'. For both of the exposure types, the zero point offset is 1-indexed and the
X (e.g., IMX) refers to the column and Y refers to the row.

MIRI MRS mode
:::::::::::::
For the MIRI MRS the SPECWCS file is in ASDF format with the following structure.

:channel: The MIRI channels in the observation, e.g. "12".
:band: The band for the observation (one of "LONG", "MEDIUM", "SHORT").
:model:
        :slice_number: The wavelength solution for each slice.
                       <slice_number> is the actual slice number (s), computed by s = channel * 100 + slice

NIRISS SOSS mode
::::::::::::::::
For NIRISS SOSS mode the SPECWCS file is in ASDF format with the following structure.

:model: A tabular model with the wavelength solution.

NIRCam Grism modes
::::::::::::::::::
For NIRCam WFSS and TSGRISM modes the SPECWCS file is in ASDF format with the following structure:

:displ: The wavelength transform models
:dispx: The x-dispersion models
:dispy: The y-dispersion models
:invdispx: The inverse x-dispersion models
:invdispy: The inverse y-dispersion models
:invdispl: The inverse wavelength transform models
:orders: a list of order numbers that the models relate to, in the same order as the models

NIRISS WFSS mode
::::::::::::::::
For NIRISS WFSS mode the SPECWCS file is in ASDF format with the following structure:

:displ: The wavelength transform models
:dispx: The x-dispersion models
:dispy: The y-dispersion models
:invdispx: The inverse x-dispersion models
:invdispl: The inverse wavelength transform models
:fwcpos_ref: The reference filter wheel position in degrees
:orders: a list of order numbers that the models relate to, in the same order as the models

