Step Arguments
==============
The ``saturation`` step has one optional argument:

``--n_pix_grow_sat`` (integer, default=1)
  The radius, in pixels, of the box to use when growing the saturation flag values
  to neighboring pixels, in order to account for charge migration (spilling). The
  box size is 2*n_pix_grow_sat+1.
