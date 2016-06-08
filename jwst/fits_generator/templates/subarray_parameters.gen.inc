/ subarray parameters

SUBARRAY  = input()    /  name of subarray used
SUBSTRT1  = int(input('COLSTART')) /  starting pixel number in the SIAS x direction
SUBSIZE1  = int(input('COLSTOP')) - int(input('COLSTART')) + 1  /  number of pixels in the SIAS x direction
SUBSTRT2  = int(input('ROWSTART'))  /  starting pixel number in the SIAS y direction
SUBSIZE2  = int(input('ROWSTOP')) - int(input('ROWSTART')) + 1   /  number of pixels in the SIAS y direction
NXLIGHT   = '#TODO'     /  number of light sensitive x values (columns)
FASTAXIS  = input()   /  Direction of fast readout (+/-1 X direction, +/-2 Y direction)
SLOWAXIS  = input()   /  Direction of slow readout (+/-1 X direction, +/-2 Y direction)
