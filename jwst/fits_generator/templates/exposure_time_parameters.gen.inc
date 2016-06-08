/ Exposure time parameters

NSAMPLES  = input('NSAMPLE') / number of A/D samples per pixel
NFRAMES   = input('NFRAME') / number of frames in group
GROUPGAP  = input() / number of frames dropped between groups
TSAMPLE   = input() / delta time between samples in microseconds
NRESET    =   1     / number of resets between integrations
NRSTSTRT  =   1     / number of extra resets at start of exposure
TFRAME    = input() / [seconds] time between frames
TGROUP    = (input('GROUPGAP') + input('NFRAME'))*input('TFRAME') / [seconds] time between groups
EFFINTTM  = ((input('NGROUP')-1)*(input('NFRAME')+input('GROUPGAP'))-input('GROUPGAP'))*input('TFRAME') / [seconds] effective integration time
EFFEXPTM  = ((input('NGROUP')-1)*(input('NFRAME')+input('GROUPGAP'))-input('GROUPGAP'))*input('TFRAME')*input('NINT') / [seconds] effective exposure time
CHRGTIME  = -1.0    / [seconds] charge accumulation time
DURATION  = -1.0    / [seconds] total duration of exposure

