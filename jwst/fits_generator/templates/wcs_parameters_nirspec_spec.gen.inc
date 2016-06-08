
/ WCS parameters

WCSAXES  = 3       /      number of World Coordinate System axes
CRPIX1   = 0.5*(input('NAXIS1')) /      x-coordinate of the reference pixel
CRPIX2   = 0.5*(input('NAXIS2')) /      y-coordinate of the reference pixel
CRVAL1   = 5.3196 /    RA at the reference pixel (degrees)
CRVAL2   = -72.986050 /    Dec at the reference pixel (degrees)
CRVAL3   = 2.5      /     Wavelength at the reference pixel (microns)
CTYPE1   = 'RA---TAN'  /  first axis coordinate type
CTYPE2   = 'DEC--TAN'  /  second axis coordinate type
CTYPE3   = 'WAVE'      /  third axis coordinate type
CUNIT1   = 'deg'       /  units for first axis
CUNIT2   = 'deg'       /  units for seconds axis
CUNIT3   = 'micron'    /  units for third axis
CDELT1   = 0.12   /      increment per pixel, axis 1
CDELT2   = 0.12   /      increment per pixel, axis 2
CDELT3   = 0.000672 /      increment per pixel, axis 3
PC1_1    = 1.0   /      linear transformation matrix element
PC1_2    = 0.0   /      linear transformation matrix element
PC2_1    = 0.0   /      linear transformation matrix element
PC2_2    = 1.0   /      linear transformation matrix element
PC3_1    = 1.0   /      linear transformation matrix element
PC3_2    = 0.0   /      linear transformation matrix element
S_REGION = 'N/A' /      spatial extent of the observation
