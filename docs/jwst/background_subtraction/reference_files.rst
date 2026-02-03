Reference Files
===============
The background image subtraction step uses reference files only when
processing Wide-Field Slitless Spectroscopy (WFSS) exposures. Two reference
files are used for WFSS mode:

* :ref:`BKG <bkg_reffile>`
* :ref:`WAVELENGTHRANGE <bg_wlrange_reffile>`

The WAVELENGTHRANGE reference file is used in the process of determining the
locations of source spectra in the image, and conversely the image areas
that contain only background signal.

.. include:: ../references_general/bkg_reffile.inc

.. _bg_wlrange_reffile:

.. include:: ../references_general/wavelengthrange_reffile.inc
