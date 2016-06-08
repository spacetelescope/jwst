/ Exposure times

EXPSTART = toMJD(input('DATE-OBS') + 'T' + input('TIME-OBS'))  / UTC exposure start time (MJD)
EXPMID   =  toMJD(input('DATE-OBS') + 'T' + input('TIME-OBS')) / UTC exposure mid time (MJD)
EXPEND   = toMJD(input('DATE-END') + 'T' + input('TIME-END')) / UTC exposure end time (MJD)

