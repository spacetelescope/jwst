/ NIRSpec configuration (NIRSpec only)

FILTER    = input('FWA_POS') or 'CLEAR'   /  name of the filter element used
GRATING   = input('GWA_POS') or 'N/A'  /     name of the grating element used
GWA_XTIL  = input('GWA_XTIL') or 0.0  /       grating x tilt
GWA_YTIL  = input('GWA_YTIL') or 0.0  /       grating y tilt
GWA_TILT  = input('GWA_TTIL') or 0.0 /    GWA temperature ave [K]
FXD_SLIT  = 'NONE'            /       name of fixed slit aperture used
MSASTATE  = 'CONFIGURED'    /  state of MSA: all open, all closed, configured
FOCUSPOS  =   0.0              /       [mm] focus position for NIRSpec

