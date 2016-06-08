/ Exposure parameters

DETECTOR = 'NRC' + input('MODULE') + str(input('SCA'))  / name of detector used to acquire data
NINTS    = input() or 1 / number of integrations within exposure
NGROUPS  = input() or input('NGROUP') / number of groups within integration
ZEROFRAM = input('ZROFRAME') / T if a zero frame was read separately
READPATT = input('READOUT') / readout pattern
DATAPROB = F                 / T if science telemetry indicated any problems with the data

