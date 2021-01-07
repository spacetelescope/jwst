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

.. include:: ../references_general/wavelengthrange_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for WAVELENGTHRANGE
++++++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in WAVELENGTHRANGE reference files

=========  ==============================
Keyword    Data Model Name
=========  ==============================
EXP_TYPE   model.meta.exposure.type
=========  ==============================

These keywords are used as CRDS selectors

.. include:: ../references_general/wavelengthrange_selection.inc


Reference File Format
+++++++++++++++++++++
WAVELENGTHRANGE reference files are in ASDF format, with the format and contents
specified by the `~jwst.datamodels.WavelengthrangeModel` data model schema.
The exact content varies a bit depending on instrument mode, as explained in the
following sections.

MIRI MRS
::::::::
For MIRI MRS, the WAVELENGTHRANGE file consists of two fields that define the
wavelength range for each combination of channel and band.

:channels: An ordered list of all possible channel and band combinations for MIRI MRS, e.g. "1SHORT".
:wavelengthrange: An ordered list of (lambda_min, lambda_max) for each item in the list above

NIRSpec
:::::::
For NIRSpec, the WAVELENGTHRANGE file is a dictionary storing information about default
wavelength range and spectral order for each combination of filter and grating.

:filter_grating:
                 :order: Default spectral order
                 :range: Default wavelength range

NIRCam WFSS, NIRCam TSGRISM, NIRISS WFSS
:::::::::::::::::::::::::::::::::::::::::
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

