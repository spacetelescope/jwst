Reference Files
===============

WCS Reference files are in the Advanced Scientific Data Format (ASDF).
The best way to create the file is to programmatically create the model and then save it to a file.
A tutorial on creating reference files in ASDF format is available at:

https://github.com/spacetelescope/jwreftools/blob/master/docs/notebooks/referece_files_asdf.ipynb

Transforms are 0-based. The forward direction is from detector to sky.

Reference file types used by ``assign_wcs``
-------------------------------------------

================================================  ====================================================  =============================
REFTYPE                                           Description                                           Instruments
================================================  ====================================================  =============================
:ref:`CAMERA <camera_reffile>`                    NIRSpec Camera model                                  NIRSpec
:ref:`COLLIMATOR <collimator_reffile>`            NIRSpec Collimator Model                              NIRSpec
:ref:`DISPERSER <disperser_reffile>`              Disperser parameters                                  NIRSpec
:ref:`DISTORTION <distortion_reffile>`            Spatial distortion model                              FGS, MIRI, NIRCam, NIRISS
:ref:`FILTEROFFSET <filteroffset_reffile>`        MIRI Imager fiter offsets                             MIRI, NIRCAM, NIRISS
:ref:`FORE <fore_reffile>`                        Transform through the NIRSpec FORE optics             NIRSpec
:ref:`FPA <fpa_reffile>`                          Transform in the NIRSpec FPA plane                    NIRSpec
:ref:`IFUFORE <ifufore_reffile>`                  Transform from the IFU slicer to the IFU entrance     NIRSpec
:ref:`IFUPOST <ifupost_reffile>`                  Transform from the IFU slicer to the back of the IFU  NIRSpec
:ref:`IFUSLICER <ifuslicer_reffile>`              IFU Slicer geometric description                      NIRSpec
:ref:`MSA <msa_reffile>`                          Transformin the NIRSpec MSA plane                     NIRSpec
:ref:`OTE <ote_reffile>`                          Transform through the Optical Telescope Element       NIRSpec
:ref:`SPECWCS <specwcs_reffile>`                  Wavelength calibration models                         MIRI, NIRCam, NIRISS
:ref:`REGIONS <regions_reffile>`                  Stores location of the regions on the detector        MIRI
:ref:`WAVELENGTHRANGE <wavelengthrange_reffile>`  Typical wavelength ranges                             MIRI, NIRCam, NIRISS, NIRSpec
================================================  ====================================================  =============================
