============
Ramp Fitting
============
The Ramp Fitting step uses two reference files: Gain and Readnoise. Both are necessary for proper computation of noise estimates.

The gain values are used to temporarily convert the pixel values from units of DN to
electrons. It is assumed that the detector gain can vary from one pixel to another, 
so gain values are stored as 2-D images.  The gain is given in units of electrons/DN.

It is assumed that the read noise can vary from pixel to pixel, so the read noise is 
also stored as a 2-D image.  The values in the reference file are assumed to be per CDS 
pair of reads, as opposed to the read noise for a single read.  The read noise is given 
in units of DN.

.. toctree::
   :maxdepth: 4

   ramp_fitting_reference_files
