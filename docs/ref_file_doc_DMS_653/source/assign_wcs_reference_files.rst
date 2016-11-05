Reference files
===============

This documents describes the reference files used by jwst.assign_wcs starting with
build 5.

In general the basic WCS keywords are in the primary header and the distortion
and spectral models are stored in reference files. All files (with one exception) are in
`ASDF <http://asdf-standard.readthedocs.org/en/latest/>`__  format.

For each observing mode, determined by the value of EXP_TYPE in the science header,
assign_wcs collects all reference files and creates a pipeline of transforms which
may include intermediate coordinate frames and the corresponding transformations
between them.

The first section lists all reference types used by assign_wcs. This corresponds to the
REFTYPE property in CRDS and the reference file. Different observing modes use different
reference types. The following sections list the reference types used by each observing mode
and how they are used.

List of reference types used by assign_wcs
------------------------------------------



===================    ==========================================================
reftype                                     description
===================    ==========================================================
**camera**             NIRSPEC Camera model
**collimator**         NIRSPEC Collimator Model
**disperser**          Disperser parameters
**distortion**         Spatial distortion model
**filteroffset**       MIRI Imager fiter offsets
**fore**               Transform through the NIRSPEC FORE optics
**fpa**                Transform in the NIRSPEC FPA plane
**msa**                Transformin the NIRSPEC MSA plane
**ote**                Transform through the Optical Telescope Element
**specwcs**            Wavelength calibration models
**regions**            Stores location of the regions on the detector
**v2v3**               Transform from MIRI instrument focal plane to V2V3 plane
**wavelengthrange**    Typical wavelength ranges
===================    ==========================================================




Observing modes supported in build 5(7?)
----------------------------------------

:MIR_IMAGE:

  | reftypes: *distortion*, *filteroffset*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, focal, CelestialFrame
  | Implements: CDP3 reference data delivery, MIRI-TN-00070-ATC_Imager_distortion_CDP_Iss5.pdf


:MIR_LRS-FIXEDSLIT, MIR_LRS-SLITLESS:

  | reftypes: *specwcs*
  | CRDS rmap rules: SUBARRAY.name: GENERIC
  | WCS pipeline coordinate frames: detector, CompositeFrame
  | Implements: CDP4 reference data delivery, MIRI-TR-10020-MPI-Calibration-Data-Description_LRSPSFDistWave_v4.0.pdf


:MIR_MRS:

  | reftypes: *distortion8, *specwcs*, *v2v3*, *wavelengthrange*
  | CRDS rmap rules: EXP_TYPE, DETECTOR, CHANNEL, BAND
  | WCS pipeline coordinate frames: detector, CompositeFrame
  | Implements: CDP4 reference data delivery, MIRI-TN-00001-ETH_Iss1-3_Calibrationproduct_MRS_d2c.pdf


:NRC_IMAGE:

  | reftypes: *distortion*
  | CRDS rmap rules: EXP_TYPE, DETECTOR, CHANNEL, BAND
  | WCS pipeline coordinate frames: detector, CelestialFrame
  | Implements: Distortion file created from TEL team data.

:NIS_IMAGE:

  | reftypes: *distortion*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, CelestialFrame
  | Implements: reference file provided by NIRISS team


:NRS_FIXEDSLIT:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, before_gwa, msa, ote, world
  | Implements: reference file provided by NIRISS team






