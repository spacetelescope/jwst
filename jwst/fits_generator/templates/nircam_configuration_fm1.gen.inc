/ NIRCam configuration (NIRCam only)

MODULE   = input('DETECTOR')[3]   /        NIRCam module: A or B
CHANNEL  = {'L': 'LONG',  \
            '1': 'SHORT', \
            '2': 'SHORT', \
            '3': 'SHORT', \
            '4': 'SHORT', \
           } [input('DETECTOR')[4]]    /       NIRCam channel: long or short
PUPIL    = input('PUPIL_P') or input('PUPIL')   /        name of the pupil element used
FILTER   = input('FILTER_P') or input('FILTER')  /        name of the filter element used
PILIN    =   F      /        pupil imaging lens in optical path: T/F

