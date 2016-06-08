/ NIRCam configuration (NIRCam only)

MODULE   = input()   /        NIRCam module: A or B
CHANNEL  = {'L': 'LONG',  \
            'S': 'SHORT', \
           } [input('DETECTOR')[0]]    /       NIRCam channel: long or short
PUPIL    = input('PUPIL')   /        name of the pupil element used
FILTER   = input('FILTER')  /        name of the filter element used
PILIN    = input('PIL_MODE')  /        pupil imaging lens in optical path: T/F

