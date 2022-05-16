Description
===========

:Class: `jwst.white_light.WhiteLightStep`
:Alias: white_light

Overview
--------
The ``white_light`` step sums the spectroscopic flux over all
wavelengths in each integration of a multi-integration extracted
spectrum product to produce an integrated ("white") flux as a
function of time for the target. This is to be applied to the ``_x1dints``
product in a spectroscopic Time-Series Observation (TSO), as part of
the :ref:`calwebb_tso3 <calwebb_tso3>` pipeline. Minimum and maximum
wavelengths may be provided to limit the summation to specified
wavelength bounds, with limits inclusive.

Input details
-------------
The input should be in the form of an ``_x1dints`` product, which contains
extracted spectra from multiple integrations for a given target.

Algorithm
---------
The algorithm performs a simple sum of the flux values over all
wavelengths for each extracted spectrum contained in the input product.
If provided, ``min_wavelength`` and ``max_wavelength`` will modify the
bounds of the sum to the specified bounds.

Output product
--------------
The output product is a table of time vs. integrated flux values, stored
in the form of a ASCII ECSV (Extended Comma-Separated Value) file.
The product type suffix is ``_whtlt``.
