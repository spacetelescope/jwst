WCS reference file information per EXP_TYPE
===========================================


:FGS_IMAGE, FGS_FOCUS, FGS_SKYFLAT, FGS_INTFLAT:

  | reftypes: *distortion*
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: reference file provided by NIRISS team

:MIR_IMAGE, MIR_TACQ, MIR_LYOT, MIR4QPM, MIR_CORONCAL:

  | reftypes: *distortion*, *filteroffset*
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: CDP6 reference data delivery, MIRI-TN-00070-ATC_Imager_distortion_CDP_Iss5.pdf

:MIR_LRS-FIXEDSLIT, MIR_LRS-SLITLESS:

  | reftypes: *specwcs*, *distortion*
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: CDP4 reference data delivery, MIRI-TR-10020-MPI-Calibration-Data-Description_LRSPSFDistWave_v4.0.pdf

:MIR_MRS:

  | reftypes: *distortion*, *specwcs*, *v2v3*, *wavelengthrange*, *regions*
  | WCS pipeline coordinate frames: detector, miri_focal, xyan, v2v3, world
  | Implements: CDP4 reference data delivery, MIRI-TN-00001-ETH_Iss1-3_Calibrationproduct_MRS_d2c.pdf

:NRC_IMAGE, NRC_TSIMAGE, NRC_FOCUS, NRC_TACONFIRM, NRC_TACQ:

  | reftypes: *distortion*, *filteroffset*
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: Distortion file created from TEL team data.

:NRC_WFSS, NRC_TSGRISM:
  | reftypes: *specwcs*, *distortion*, *filteroffset*
  | WCS pipeline coordinate frames: grism_detector, detector, v2v3, world
  | Implements: reference files provided by NIRCam team

:NIS_IMAGE, NIS_TACQ, NIS_TACONFIRM, NIS_FOCUS:

  | reftypes: *distortion*, *filteroffset*
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: reference file provided by NIRISS team

:NIS_WFSS:
  | reftypes: *specwcs*, *distortion*, *filteroffset*
  | WCS pipeline coordinate frames: grism_detector, detector, v2v3, world
  | Implements: reference files provided by NIRISS team

:NIS_SOSS:

  | reftypes: *distortion*, *specwcs*
  | WCS pipeline coordinate frames: detector, v2v3, world
  | Implements: reference files provided by NIRISS team

:NRS_FIXEDSLIT, NRS_MSASPEC, NRS_LAMP, NRS_BRIGHTOBJ:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*
  | WCS pipeline coordinate frames: detector, sca, bgwa, slit_frame, msa_frame, ote, v2v3, world
  | Implements: CDP 3 delivery

:NRS_IFU:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*,
  | *ifufore*, *ifuslicer*, *ifupost*
  | WCS pipeline coordinate frames: detector, sca, bgwa, slit_frame, msa_frame, ote, v2v3, world
  | Implements: CDP 3 delivery

:NRS_IMAGING, NRS_MIMF, NRS_BOTA, NRS_CONFIRM, NRS_TACONFIRM, NRS_TASLIT, NRS_TACQ:

  | reftypes: *fpa*, *camera*, *disperser*, *collimator*, *msa*, *wavelengthrange*, *fore*, *ote*
  | WCS pipeline coordinate frames: detector, sca, bgwa, slit_frame, msa_frame, ote, v2v3, world
  | Implements: CDP 3 delivery
