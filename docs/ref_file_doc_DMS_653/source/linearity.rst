=========
Linearity
=========


The linearity correction corrects the integrated counts in the science images for the non-linear response of the detector.
The correction is applied pixel-by-pixel, group-by-group, integration-by-integration within a science exposure.
The correction is represented by an nth-order polynomial for each pixel in the detector, with n+1 arrays of coefficients 
read from the linearity reference file.  The values from the linearity reference file DQ array are propagated into the
PIXELDQ array of the input science exposure using a bitwise OR operation.


.. toctree::
   :maxdepth: 4

   linearity_reference_files
