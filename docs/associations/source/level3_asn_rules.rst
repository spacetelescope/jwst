Level3 Associations: Rules
``````````````````````````
.. _asn-association-types:

Association Types
=================

Each Level3 association is intended to make a specific science
product. The type of science product is indicated by the `ATYPE` field
in the association file name (see :ref:`asn-DMS-naming`), and in the `asn_type` meta
keyword of the association itself (see :ref:`asn-association-meta-keywords`).

The pipeline uses this type as the key to indicate which Level 3
pipeline module to use to process this association.

The current association types are:

  * `image`: suitable for CALIMAGE3 processing
  * `spec`: suitable for CALSPECE3 processing
  * `wfs`: Wave front sensing data, used by `wfs_combine`
  * `ami`: Aperture Mask Interferometry
  * `coron`: Coronography
  * `tso`: Time-series Observations
  * `wfss`: Wide-Field Slitless Spectroscopy
  * `nrsifu`: NIRSpec IFU Spectroscopy
    
    This is different from just spectroscopy because of how NIRSpec
    is arranged. The wavelength/dispersion is determined by the
    Prism/Grating wheel. The Filter wheel does not modified
    wavelenght/dispersion, thus allowing exposures in multiple
    wavelengths to be combined.


Rules
=====
