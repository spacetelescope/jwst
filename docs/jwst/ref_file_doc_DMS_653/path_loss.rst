=========
Path Loss
=========

The pathloss reference file gives the path loss correction as a function of wavelength.  
There are two types of pathloss calibrations performed: for point sources and for uniform sources.

The point source entry in the reference file is a 3-d array with the pathloss correction as a 
function of wavelength and decenter within the aperture.  The pathloss correction interpolates 
the 3-d array at the location of a point source to provide a 1-d array of pathloss vs. wavelength.  
This 1-d array is attached to the data model in the pathloss_pointsource attribute, with corresponding 
wavelength array in the wavelength_pointsource attribute.


The uniform source entry has a 1-d array of pathloss vs. wavelength.  This array is attached 
to the data model in the pathloss_uniformsource attribute, along with the wavelength array in 
the wavelength_uniformsource attribute.


.. toctree::
   :maxdepth: 4

   path_loss_reference_files
