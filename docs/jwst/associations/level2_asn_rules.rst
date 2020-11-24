.. asn-level2-rules:

Level 2 Associations: Rules
===========================

The following table describes exactly which exposures will go
through any type of Level 2b processing, such as Spec2Pipeline or
Image2Pipeline.

.. list-table:: Exposure Modes for Level 2 Processing
   :widths: 20 20 20 20
   :header-rows: 1

   * - EXP_TYPE
     - Association Exposure Type
     - Specials
     - Level2b?
   * - NRC_IMAGE
     - science
     - n/a
     - YES
   * - NRC_GRISM/NRC_WFSS
     - science
     - n/a
     - YES
   * - NRC_TACQ
     - target_acquisition
     - n/a
     - YES
   * - NRC_TACONFIRM
     - target_acquisition
     - n/a
     - YES
   * - NRC_CORON
     - science
     - n/a
     - YES
   * - NRC_CORON
     - psf
     - PSF
     - YES
   * - NRC_TSIMAGE
     - science
     - n/a
     - YES
   * - NRC_TSGRISM
     - science
     - n/a
     - YES
   * - NRC_FOCUS
     - science
     - n/a/
     - YES
   * - NRC_DARK
     - dark
     - n/a
     - NO
   * - NRC_FLAT
     - flat
     - n/a
     - NO
   * - NRC_LED
     - science
     - n/a
     - NO
   * -
     -
     -
     -
   * - MIR_IMAGE
     - science
     - n/a
     - YES
   * - MIR_TACQ
     - target_acquisition
     - n/a
     - YES
   * - MIR_LYOT
     - science
     - n/a
     - YES
   * - MIR_LYOT
     - psf
     - PSF
     - YES
   * - MIR_4QPM
     - science
     - n/a
     - YES
   * - MIR_4QPM
     - psf
     - PSF
     - YES
   * - MIR_LRS-FIXEDSLIT
     - science
     - n/a
     - YES
   * - MIR_LRS-FIXEDSLIT
     - background
     - BACKGROUND
     - YES
   * - MIR_LRS-SLITLESS
     - science
     - n/a
     - YES
   * - MIR_LRS-SLITLESS
     - background
     - BACKGROUND
     - YES
   * - MIR_MRS
     - science
     - n/a
     - YES
   * - MIR_MRS
     - background
     - BACKGROUND
     - YES
   * - MIR_DARKIMG
     - dark
     - n/a
     - NO
   * - MIR_DARKMRS
     - dark
     - n/a
     - NO
   * - MIR_FLATIMAGE
     - flat
     - n/a
     - NO
   * - MIR_FLATIMAGE-EXT
     - flat
     - n/a
     - NO
   * - MIR_FLATMRS
     - flat
     - n/a
     - NO
   * - MIR_FLATMRS-EXT
     - flat
     - n/a
     - NO
   * - MIR_CORONCAL
     - science
     - n/a
     - YES
   * -
     -
     -
     - 
   * - NRS_WATA
     - target_acquisition
     - n/a
     - YES
   * - NRS_MSATA
     - target_acquisition
     - n/a
     - YES
   * - NRS_TACONFIRM
     - target_acquisition
     - n/a
     - YES
   * - NRS_CONFIRM
     - science
     - n/a
     - YES
   * - NRS_FIXEDSLIT
     - science
     - n/a
     - YES
   * - NRS_FIXEDSLIT
     - background
     - BACKGROUND
     - YES
   * - NRS_AUTOWAVE
     - nrs_autowave
     - n/a
     - YES
   * - NRS_IFU
     - science
     - n/a
     - YES
   * - NRS_IFU
     - imprint
     - IMPRINT
     - YES
   * - NRS_IFU
     - background
     - BACKGROUND
     - YES
   * - NRS_IMAGE
     - science
     - n/a
     - YES
   * - NRS_MSASPEC
     - science
     - n/a
     - YES
   * - NRS_MSASPEC
     - imprint
     - IMPRINT
     - YES
   * - NRS_AUTOFLAT
     - nrs_autoflat
     - n/a
     - YES
   * - NRS_FOCUS
     - science
     - n/a
     - YES
   * - NRS_DARK
     - dark
     - n/a
     - NO
   * - NRS_LAMP
     - science
     - n/a
     - YES<sup>1</sup>
   * - NRS_BRIGHTOBJ
     - science
     - n/a
     - YES
   * - NRS_MIMF
     - science
     - n/a
     - YES
   * - NRS_VERIFY
     - science
     - n/a
     - YES<sup>1</sup>
   * -
     -
     -
     - 
   * - NIS_IMAGE
     - science
     - n/a
     - YES
   * - NIS_WFSS
     - science
     - n/a
     - YES
   * - NIS_TACQ
     - target_acquisition
     - n/a
     - YES
   * - NIS_TACONFIRM
     - target_acquisition
     - n/a
     - YES
   * - NIS_SOSS
     - science
     - n/a
     - YES
   * - NIS_AMI
     - science
     - n/a
     - YES
   * - NIS_AMI
     - psf
     - PSF
     - YES
   * - NIS_FOCUS
     - science
     - n/a
     - YES
   * - NIS_DARK
     - science
     - n/a
     - NO
   * - NIS_LAMP
     - science
     - n/a
     - NO
   * - NIS_EXTCAL
     - science
     - n/a
     - NO
   * -
     -
     -
     - 
   * - FGS_IMAGE
     - science
     - n/a
     - YES
   * - FGS_FOCUS
     - science
     - n/a
     - YES
   * - FGS_SKYFLAT
     - flat
     - n/a
     - NO
   * - FGS_INTFLAT
     - flat
     - n/a
     - NO
   * - FGS_DARK
     - dark
     - n/a
     - NO
   * - FGS_ID-STACK
     - tracking
     - n/a
     - NO
   * - FGS_ID-IMAGE
     - tracking
     - n/a
     - NO
   * - FGS_ACQ1
     - tracking
     - n/a
     - NO
   * - FGS_ACQ2
     - tracking
     - n/a
     - NO
   * - FGS_TRACK
     - tracking
     - n/a
     - NO
   * - FGS_FINEGUIDE
     - tracking
     - n/a
     - NO

History
-------

The original content of this page is from `github issue #1188`_.

.. _github issue #1188: https://github.com/spacetelescope/jwst/issues/1188
