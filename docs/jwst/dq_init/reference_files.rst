Reference File Types
====================
The Data Quality Initialization step uses a MASK reference file.

.. include:: mask_selection.rst

MASK Reference File Format
--------------------------
The MASK reference file is a FITS file with a primary HDU, 1 IMAGE extension
HDU and 1 BINTABLE extension.

The primary data array is assumed to be empty.  The MASK data are stored in the
first IMAGE extension, which shall have EXTNAME='DQ'. The data array in this
extension has integer data type and is 2-D, with dimensions equal to the number
of columns and rows in a full frame raw readout for the given detector,
including reference pixels. Note that this does not include the reference
output for MIRI detectors.

The BINTABLE extension contains the bit assignments used in the DQ array.

.. include:: ../includes/dq_def.rst

.. include:: fgs_mask.rst
.. include:: miri_mask.rst
.. include:: nircam_mask.rst
.. include:: niriss_mask.rst
.. include:: nirspec_mask.rst

