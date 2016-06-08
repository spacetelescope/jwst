/ NIRISS configuration (NIRISS only)

FILTER   =   input('FWCCRFIL')   /   name of the filter element used
PUPIL    =   {'CLEAR': 'CLEARP', \
              'F090W': 'F090W',  \
              'F115W': 'F115W',  \
              'F150W': 'F150W',  \
              'F200W': 'F200W',  \
              'UNKNOWN': 'UNKNOWN',  \
              } [input('PWCCRPUP')]   /   name of the pupil element used
FOCUSPOS =     0.0     /   [mm] focus position for NIRISS

