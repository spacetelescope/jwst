SIMPLE   = T       / Data conform to FITS standard
BITPIX   = 0       / Bits per data value
NAXIS    = 0       / Number of data array dimensions
EXTEND   = T       / File may contain standard extensions
NEXTEND  = 3       / Number of standard extensions

TELESCOP = 'JWST'  / Telescope used to acquire data
INSTRUME = input() / Identifier for instrument used to acquire data
RADESYS  = 'ICRS'  / Coordinate reference frame for RA and Dec

DATE     = now()   / UTC date file created
ORIGIN   = 'STScI' / Institution responsible for creating FITS file
FILENAME = '#TODO' / Name of file
FILETYPE = 'UNCALIBRATED' / Type of data found in file
DPSW_VER = version() / Data processing software version number

