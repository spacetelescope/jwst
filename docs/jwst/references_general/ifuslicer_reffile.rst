:orphan:

.. _ifuslicer_reffile:
  
IFUSLICER Reference File (NIRSpec only)
---------------------------------------

:REFTYPE: IFUSLICER
:Data model: `~jwst.datamodels.IFUSlicerModel`

Reference Selection Keywords for IFUSLICER
++++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate IFUSLICER references based on the following keywords.
IFUSLICER is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ======================================
Instrument Keywords
========== ======================================
NIRSpec    INSTRUME, EXP_TYPE, DATE-OBS, TIME-OBS
========== ======================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The IFUSLICER stores information about the metrology of the IFU slicer - relative
positioning and size of the aperture of each individual slicer and the absolute
reference with respect to the center of the field of view.
The reference file contains two fields - “data” and “model”.
The “data” field is an array with 30 rows pertaining to the 30 slices and the columns are

:data: Array with reference data for each slicer. It has 5 columns

          NO
            Slice number (0 - 29)
          x_center
            X coordinate of the center (in meters)
          y_center
            Y coordinate of the center (in meters)
          x_size
            X size of teh aperture (in meters)
          y_size
            Y size of the aperture (in meters)
:model: Transform from relative positions within the IFU slicer to absolute positions
        within the field of view. It's a combination of shifts and rotation.
