.. asn-level2-rules:

Stage 2 Associations: Rules
===========================

The following table describes exactly which exposures will go
through any type of stage 2 processing, such as Spec2Pipeline or
Image2Pipeline.

.. list-table:: Exposure Modes for Stage 2 Processing
   :widths: 20 20 20 20
   :header-rows: 1

   * - EXP_TYPE
     - Member Exposure Type
     - Specials
     - Association Type
   * - FGS_ACQ1
     - tracking
     - N/A
     - N/A
   * - FGS_ACQ2
     - tracking
     - N/A
     - N/A
   * - FGS_DARK
     - dark
     - N/A
     - N/A
   * - FGS_FINEGUIDE
     - tracking
     - N/A
     - N/A
   * - FGS_FOCUS
     - science
     - N/A
     - image2
   * - FGS_ID-IMAGE
     - tracking
     - N/A
     - N/A
   * - FGS_ID-STACK
     - tracking
     - N/A
     - N/A
   * - FGS_IMAGE
     - science
     - N/A
     - image2
   * - FGS_INTFLAT
     - flat
     - N/A
     - N/A
   * - FGS_SKYFLAT
     - flat
     - N/A
     - N/A
   * - FGS_TRACK
     - tracking
     - N/A
     - N/A
   * -
     -
     -
     - 
   * - MIR_4QPM
     - psf
     - PSF
     - image2
   * - MIR_4QPM
     - science
     - N/A
     - image2
   * - MIR_CORONCAL
     - science
     - N/A
     - image2
   * - MIR_DARKIMG
     - dark
     - N/A
     - N/A
   * - MIR_DARKMRS
     - dark
     - N/A
     - N/A
   * - MIR_FLATIMAGE
     - flat
     - N/A
     - N/A
   * - MIR_FLATIMAGE-EXT
     - flat
     - N/A
     - N/A
   * - MIR_FLATMRS
     - flat
     - N/A
     - N/A
   * - MIR_FLATMRS-EXT
     - flat
     - N/A
     - N/A
   * - MIR_IMAGE
     - science
     - N/A
     - image2
   * - MIR_LRS-FIXEDSLIT
     - background
     - BACKGROUND
     - spec2
   * - MIR_LRS-FIXEDSLIT
     - science
     - N/A
     - spec2
   * - MIR_LRS-SLITLESS
     - background
     - BACKGROUND
     - spec2
   * - MIR_LRS-SLITLESS
     - science
     - N/A
     - spec2
   * - MIR_LYOT
     - psf
     - PSF
     - image2
   * - MIR_LYOT
     - science
     - N/A
     - image2
   * - MIR_MRS
     - background
     - BACKGROUND
     - spec2
   * - MIR_MRS
     - science
     - N/A
     - spec2
   * - MIR_TACQ
     - target_acquisition
     - N/A
     - image2
   * -
     -
     -
     - 
   * - NIS_AMI
     - psf
     - PSF
     - image2
   * - NIS_AMI
     - science
     - N/A
     - image2
   * - NIS_DARK
     - science
     - N/A
     - N/A
   * - NIS_EXTCAL
     - science
     - N/A
     - N/A
   * - NIS_FOCUS
     - science
     - N/A
     - image2
   * - NIS_IMAGE
     - science
     - N/A
     - images
   * - NIS_LAMP
     - science
     - N/A
     - N/A
   * - NIS_SOSS
     - science
     - N/A
     - spec2
   * - NIS_TACONFIRM
     - target_acquisition
     - N/A
     - image2
   * - NIS_TACQ
     - target_acquisition
     - N/A
     - image2
   * - NIS_WFSS
     - science
     - N/A
     - spec2
   * -
     -
     -
     - 
   * - NRC_CORON
     - psf
     - PSF
     - image2
   * - NRC_CORON
     - science
     - N/A
     - image2
   * - NRC_DARK
     - dark
     - N/A
     - N/A
   * - NRC_FLAT
     - flat
     - N/A
     - N/A
   * - NRC_FOCUS
     - science
     - N/A/
     - image2
   * - NRC_GRISM
     - science
     - N/A
     - N/A
   * - NRC_IMAGE
     - science
     - N/A
     - image2
   * - NRC_LED
     - science
     - N/A
     - N/A
   * - NRC_TACONFIRM
     - target_acquisition
     - N/A
     - image2
   * - NRC_TACQ
     - target_acquisition
     - N/A
     - image2
   * - NRC_TSGRISM
     - science
     - N/A
     - tso-spec2
   * - NRC_TSIMAGE
     - science
     - N/A
     - tso-image2
   * - NRC_WFSS
     - science
     - N/A
     - spec2
   * -
     -
     -
     -
   * - NRS_AUTOFLAT
     - nrs_autoflat
     - N/A
     - image2
   * - NRS_AUTOWAVE
     - nrs_autowave
     - N/A
     - image2
   * - NRS_BRIGHTOBJ
     - science
     - N/A
     - spec2
   * - NRS_CONFIRM
     - science
     - N/A
     - image2
   * - NRS_DARK
     - dark
     - N/A
     - N/A
   * - NRS_FIXEDSLIT
     - background
     - BACKGROUND
     - spec2
   * - NRS_FIXEDSLIT
     - science
     - N/A
     - spec2
   * - NRS_FOCUS
     - science
     - N/A
     - image2
   * - NRS_IFU
     - background
     - BACKGROUND
     - spec2
   * - NRS_IFU
     - imprint
     - IMPRINT
     - spec2
   * - NRS_IFU
     - science
     - N/A
     - spec2
   * - NRS_IMAGE
     - science
     - N/A
     - image2
   * - NRS_LAMP [#f1]_
     - science
     - N/A
     - nrslamp-spec2
   * - NRS_MIMF
     - science
     - N/A
     - wfs-image2
   * - NRS_MSASPEC
     - imprint
     - IMPRINT
     - spec2
   * - NRS_MSASPEC
     - science
     - N/A
     - spec2
   * - NRS_MSATA
     - target_acquisition
     - N/A
     - image2
   * - NRS_TACONFIRM
     - target_acquisition
     - N/A
     - image2
   * - NRS_VERIFY
     - science
     - N/A
     - image2
   * - NRS_WATA
     - target_acquisition
     - N/A
     - image2

Footnotes
---------

.. [#f1] Association creation is heavily dependent upon other parameters such as ``LAMP``, ``OPMODE``, and ``GRATING``.

Notes
-----

Column definitions

- EXP_TYPE : The exposure type.
- Member Exposure Type: How the association generator will classify the exposure.
- Specials : The association rule modifications to handle the exposure.
- Association Type : :ref:`Association type <asn-jwst-association-types>` created.

More about Specials: Many exposures that are not directly science, such as
backgrounds, are primarily used as auxiliary members for other science products.
However, they are also often calibrated as if they were science products
themselves. In these situations, a special association rule is created to
produce the necessary associations.

History
-------

The original content of this page is from `github issue #1188`_.

.. _github issue #1188: https://github.com/spacetelescope/jwst/issues/1188
