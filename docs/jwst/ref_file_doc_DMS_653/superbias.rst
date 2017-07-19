=========
Superbias
=========
The superbias subtraction step removes the fixed detector bias from a
science data set by subtracting a superbias reference image.
The 2-D superbias reference image is subtracted from every group in every
integration of the input science ramp data. Any NaN's that might be present
in the superbias image are set to a value of zero before being subtracted
from the science data, such that those pixels effectively receive no correction.
The DQ array from the superbias reference file is combined with the science
exposure PIXELDQ array using a bit-wise OR operation.

.. toctree::
   :maxdepth: 4

   superbias_reference_files
