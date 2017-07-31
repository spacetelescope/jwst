===========
Extract_2d
===========


The extract_2d step extracts a 2-D cutout for each spectrum in an exposure.
It is saved as a "SCI" extension. It also computes and saves the wavelengths
in a separate extension with EXTNAME "WAVELENGTH". It works on Nirspec MSA and
fixed slits, as well as on NIRISS and NIRCAM slitless observations. Point source
Nirspec wavelengths are (optionally) corrected for 
An optional wavelength zero-point correction is applied to Nirspec point source
observations when the source is not centered in the slit. The data for the correction
is saved in a WAVECORR reference file.

.. toctree::
   :maxdepth: 4

   extract_2d_reference_files
