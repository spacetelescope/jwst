Step Arguments
==============
The ``saturation`` step has one optional argument:

``--n_pix_grow_sat`` (integer, default=1)
  The distance to use when growing saturation flag values to neighboring pixels,
  in order to account for charge migration (spilling). The total region size is
  2*n_pix_grow_sat+1 pixels, centered on the primary pixel.
