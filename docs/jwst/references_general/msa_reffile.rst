:orphan:

.. _msa_reffile:
  
MSA Reference File (NIRSpec only)
---------------------------------

:REFTYPE: MSA
:Data model: `~jwst.datamodels.MSAModel`

Reference Selection Keywords for MSA
++++++++++++++++++++++++++++++++++++
CRDS selects appropriate MSA references based on the following keywords.
MSA is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ======================================
Instrument Keywords
========== ======================================
NIRSpec    INSTRUME, EXP_TYPE, DATE-OBS, TIME-OBS
========== ======================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The MSA reference file contains information on the metrology of the microshutter array and the associated fixed slits - relative positioning of each individual shutter (assumed to be rectangular)
And the absolute position of each quadrant within the MSA.

The MSA reference file has 5 fields, named

:1:
   :data: Array with reference data for each shutter in Quadrant 1.
          It has 5 columns

          NO
            Shutter number (1- 62415)
          x_center
            X coordinate of the center (in meters)
          y_center
            Y coordinate of the center (in meters)
          x_size
            X size of teh aperture (in meters)
          y_size
            Y size of the aperture (in meters)
   :model: Transform from relative positions within Quadrant 1 to absolute positions within the MSA
:2:
   :data: Array with reference data for shutters in Quadrant 2, same as in 1 above
   :model: Transform from relative positions within Quadrant 2 to absolute positions within the MSA
:3:
   :data: Array with reference data for shutters in Quadrant 3, same as in 1 above
   :model: Transform from relative positions within Quadrant 3 to absolute positions within the MSA
:4:
   :data: Array with reference data for shutters in Quadrant 4, same as in 1 above
   :model: Transform from relative positions within Quadrant 4 to absolute positions within the MSA
:5:
   :data: Reference data for the fixed slits and the IFU, same as in 1, except NO is 6 rows (1-6)
          and the mapping is 1 - S200A1, 2 - S200A1, 3 - S400A1, 4 - S200B1, 5 - S1600A1, 6 - IFU
   :model: Transform from relative positions within eac aperture to absolute positions within the MSA

