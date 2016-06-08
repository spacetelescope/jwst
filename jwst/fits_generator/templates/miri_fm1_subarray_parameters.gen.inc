/ subarray parameters

SUBARRAY  = input()    /  subarray used
SUBSTRT1  = int(input('COLSTART')) /  starting pixel number in the SIAS x direction
SUBSIZE1  = int(input('NAXIS1')) /  number of pixels in the SIAS x direction
SUBSTRT2  = int(input('ROWSTART'))  /  starting pixel number in the SIAS y direction
SUBSIZE2  = int(input('ROWSTOP')) - int(input('ROWSTART')) + 1   /  number of pixels in the SIAS y direction
LIGHTCOL  = '#TODO  '  /  number of light sensitive columns
FASTAXIS  = input()    /  Direction of fast readout
SLOWAXIS  = input()    /  Direction of slow readout
