======
Photom
======

The photom step copies flux conversion information from the photometric reference table into 
the science product.  The step searches the reference table for the row that matches the parameters of the exposure;
the row contains a scalar conversion constant, as well as optional arrays
of wavelength and relative response (as a function of wavelength).  The scalar conversion 
constant is copied into the keyword PHOTMJSR in the primary header of the science product, and, if 
the wavelength and relative response arrays are populated in the selected row, those arrays are 
copied to a table extension called “RELSENS”.

If the science data are from an imaging mode, the data from the pixel area map reference 
file will also be copied into the science data product. The 2-D data array from
the pixel area map will be copied into an image extension called “AREA”, and the 
values of the PIXAR_SR and PIXAR_A2 keywords in the photom reference table will also be
copied into keywords of the same name in the primary header.

.. toctree::
   :maxdepth: 4

   photom_reference_files
