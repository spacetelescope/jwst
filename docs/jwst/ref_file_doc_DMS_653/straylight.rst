==========
Straylight
==========

The stray-light step is applied to MIRI MRS data only, in which case a stray-light MASK reference
file is used to designate which pixels are science pixels and which pixels fall in-between the slices.
Each illuminated pixel on the array has a signal that is the sum of direct
illumination and the scattering from neighboring areas. Only the pixels located between the slices are areas of indirect
illumination. The illumination on the inter-slice pixels are used to determine a stray-light component to subtract from
each science pixel.

.. toctree::
   :maxdepth: 4

   straylight_reference_files
