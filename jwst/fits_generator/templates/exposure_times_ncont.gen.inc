/ Exposure times

EXPSTART = toMJD(input('DATE')[6:10] + '-' + input('DATE')[0:2] + '-' + input('DATE')[3:5] + 'T' + input('STARTEXP'))  / UTC exposure start time (MJD)
EXPMID   =  toMJD(input('DATE')[6:10] + '-' + input('DATE')[0:2] + '-' + input('DATE')[3:5] + 'T' + input('STARTEXP')) / UTC exposure mid time (MJD)
EXPEND   = toMJD(input('DATE')[6:10] + '-' + input('DATE')[0:2] + '-' + input('DATE')[3:5] + 'T' + input('END_EXP')) / UTC exposure end time (MJD)

