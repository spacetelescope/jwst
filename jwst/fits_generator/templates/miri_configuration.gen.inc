/ MIRI configuration (MIRI only)

FILTER    = input('FWA_POS')   /       name of the filter element used
BAND      = input('DGAA_POS')  /       spectral band setting (SHORT, MEDIUM or LONG)
CHANNEL   = {493: 'N/A', \
             494: '34', \
             495: '12'}[input('SCA_ID')]
CCCSTATE  = input('CCC_POS') / contamination control cover state

