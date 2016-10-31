Reference files and WCS keywords
================================

This documents describes the reference files used by jwst.assign_wcs in
build 7.

In general the basic WCS keywords are in the primary header and the distortion
and spectral models are stored in reference files. All files (with one exception) are in
`ASDF <http://asdf-standard.readthedocs.org/en/latest/>`__  format.

For each observing mode, determined by the value of EXP_TYPE in the science header,
assign_wcs retrieves reference files from CRDS and creates a pipeline of transforms from
input frame `detector` to a frame `v2v3`. This part of the WCS pipeline may include
intermediate coordinate frames. The basic WCS keywords are used to create
the transforms from frame `v2v3` to frame `world`.

Basic WCS keywords and the transform from `v2v3` to `world`
-----------------------------------------------------------

All JWST instruments use the following FITS header keywords to define the transform from `v2v3` to `world`:

`RA_REF`, `DEC_REF` - a fiducial point on the sky, ICRS, [deg]
`V2_REF`, `V3_REF` - a point in the V2V3 system which maps to `RA_REF`, `DEC_REF`, [arcsec]
`ROLL_REF` - local roll angle associated with each aperture, [deg]

These quantities are used to create a 3D Euler angle rotation between the V2V3 spherical system,
associated with the telescope, and a standard celestial system.

Reference files used in build 7
-------------------------------

The first section lists all reference types used by assign_wcs. This corresponds to the
REFTYPE property in CRDS and the reference file. Different observing modes use different
reference types. 
The following sections list the reference types used by each observing mode
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
**ifufore**            Transform from the IFU slicer to the IFU entrance
**ifupost**            Transform from the IFU slicer to the back of the IFU
**ifuslicer**          FU Slicer geometric description
**msa**                Transformin the NIRSPEC MSA plane
**ote**                Transform through the Optical Telescope Element
**specwcs**            Wavelength calibration models
**regions**            Stores location of the regions on the detector
**v2v3**               Transform from MIRI instrument focal plane to V2V3 plane
**wavelengthrange**    Typical wavelength ranges
===================    ==========================================================




Observing modes supported in build 7
------------------------------------

:FGS_IMAGE:

  | reftypes: *distortion*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: reference file provided by NIRISS team

:MIR_IMAGE:

  | reftypes: *distortion*, *filteroffset*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: CDP6 reference data delivery, MIRI-TN-00070-ATC_Imager_distortion_CDP_Iss5.pdf


:MIR_LRS-FIXEDSLIT, MIR_LRS-SLITLESS:

  | reftypes: *specwcs*, *distortion*
  | CRDS rmap rules: SUBARRAY.name: GENERIC
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: CDP6 reference data delivery, MIRI-TR-10020-MPI-Calibration-Data-Description_LRSPSFDistWave_v4.0.pdf


:MIR_MRS:

  | reftypes: *distortion*, *specwcs*, *v2v3*, *wavelengthrange*, *regions*
  | CRDS rmap rules: EXP_TYPE, DETECTOR, CHANNEL, BAND
  | WCS pipeline coordinate frames: detector, miri_focal, xyan, v2v3, world
  | Implements: CDP4 reference data delivery, MIRI-TN-00001-ETH_Iss1-3_Calibrationproduct_MRS_d2c.pdf

:NRC_IMAGE:

  | reftypes: *distortion*
  | CRDS rmap rules: EXP_TYPE, DETECTOR, CHANNEL, BAND
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: Distortion file created from TEL team data.

:NIS_IMAGE:

  | reftypes: *distortion*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: reference file provided by NIRISS team

:NIS_SOSS:

  | reftypes: *distortion*, *specwcs*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: reference files provided by NIRISS team

:NRS_FIXEDSLIT:
:NRS_MSASPEC:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, sca, bgwa, slit_frame, msa_frame, ote, v2v3, world
  | Implements: CDP 2 delivery

:NRS_IFU:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*,
  | *ifufore*, *ifuslicer*, *ifupost*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, sca, bgwa, slit_frame, msa_frame, ote, v2v3, world
  | Implements: CDP 2 delivery

:NRS_IMAGING:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*
  | CRDS rmap rules: EXP_TYPE
  | WCS pipeline coordinate frames: detector, sca, bgwa, slit_frame, msa_frame, ote, v2v3, world
  | Implements: CDP 2 delivery




