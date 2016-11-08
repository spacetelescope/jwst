======
Refpix
======

The refpix step corrects for bias drift.  The reference-pixel reference file contains 
frequency-dependent weights that are used to compute (in Fourier space) 
the filtered reference pixels and reference output for the reference-pixel correction scheme that 
is applied to NIRSpec data when exposures use the IRS2 readout pattern.  Only the NIRSpec IRS2 
readout format requires a reference file; no other instruments or exposure modes require a reference
file for this step.

For each sector, the correction is applied as follows:  data * alpha[i] + reference_output * beta[i].
`Alpha` and `beta` are 2-D arrays of values read from the reference file.  The first axis is 
the sector number (but only for the normal pixel data and reference pixels, not the reference output).  
The second axis has length 2048 * 712, corresponding to the time-ordered arrangement of the 
data.  `Data` is the science data for the current integration.  The shape is expected to 
be (ngroups, ny, 3200), where ngroups is the number of groups, and ny is the pixel height 
of the image.  The width 3200 of the image includes the "normal" pixel data, plus the embedded 
reference pixels, and the reference output.  `Reference_output` is the length of the reference 
output section.


.. toctree::
   :maxdepth: 4

   refpix_reference_files
