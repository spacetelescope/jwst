Description
===========

Overview
--------

The ``white_light`` step sums the spectroscopic flux over all
wavelengths in each integration of a multi-integration extracted
spectrum product to produce an integrated ("white") flux as a
function of time for a target. This is usually applied to the _x1dints
product in a spectroscopic Time-Series Observation (TSO), as part of
the ``caltso3`` pipeline.

Input details
-------------
The input should be in the form of an _x1dints product, which contains
extracted spectra from multiple integrations for a given target.

Algorithm
---------
The algorithm performs a simple sum of the flux values over all
wavelengths for each extracted spectrum contained in the input product.

Output product
--------------
The output product will be a table of time vs. integrated flux, stored
in the form of a ECSV (Extended Comma-Separated Value) file.
