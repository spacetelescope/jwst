Reference Files
===============

The ``extract_2d`` step uses WAVECORR and WAVELENGTHRANGE reference files.
The WAVECORR reference file is only used for NIRSpec fixed-slit and MOS
exposures.
The WAVELENGTHRANGE reference file is only used for NIRCam and NIRISS
Wide-Field Slitless Spectroscopy (WFSS) exposures.

WAVECORR Reference File
-----------------------

:REFTYPE: WAVECORR
:Data model: `~jwst.datamodels.WaveCorrModel`

The WAVECORR reference file contains pixel offset values as a function of
wavelength and source offset within a NIRSpec slit.
It is used when applying the NIRSpec wavelength zero-point correction to
fixed-slit (EXP_TYPE="NRS_FIXEDSLIT"), bright object TSO
(EXP_TYPE="NRS_BRIGHTOBJ"), and MSA/MOS spectra (EXP_TYPE="NRS_MSASPEC").
This is an optional correction that is turned on by default.
It can be turned off by specifying ``apply_wavecorr=False`` when running the step.

.. include:: wavecorr_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for WAVECORR
+++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in WAVECORR reference files,
because they are used as CRDS selectors
(see :ref:`wavecorr_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
+++++++++++++++++++++
WAVECORR reference files are in ASDF format, with the format and contents
specified by the `~jwst.datamodels.WaveCorrModel` data model schema.

WAVELENGTHRANGE Reference File
------------------------------

:REFTYPE: WAVELENGTHRANGE
:Data model: `~jwst.datamodels.WavelengthrangeModel`

The WAVELENGTHRANGE reference file contains information on the minimum and
maximum wavelengths of various spectroscopic modes, which can be
order-dependent. The reference data are used to construct bounding boxes
around the spectral traces produced by each object in the NIRCam and NIRISS
WFSS modes.
If a list of ``GrismObject`` is supplied, then no reference file is neccessary.

.. include:: wavelengthrange_selection.rst

.. include:: ../includes/standard_keywords.rst

Type Specific Keywords for WAVELENGTHRANGE
++++++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in WAVELENGTHRANGE reference files,
because they are used as CRDS selectors
(see :ref:`wavelengthrange_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

Reference File Format
+++++++++++++++++++++
WAVELENGTHRANGE reference files are in ASDF format, with the format and contents
specified by the `~jwst.datamodels.WavelengthrangeModel` data model schema.
For NIRCam WFSS and TSGRIM modes, as well as NIRISS WFSS mode, the WAVELENGTHRANGE
reference file contains the wavelength limits to use when calculating the minimum
and maximum dispersion extents on the detector. It also contains the default list
of orders that should be extracted for each filter.
To be consistent with other modes, and for convenience, it also lists the orders
and filters that are valid with the file.

:order: A list of orders this file covers
:wavelengthrange: A list containing the list of [order, filter, wavelength min, wavelength max]
:waverange_selector: The list of FILTER names available
:extract_orders: A list containing the list of orders to extract for each filter

