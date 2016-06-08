/ Exposure parameters

DETECTOR = input() / name of detector used to acquire data
NINTS    = int(input('NINT')) / number of integrations within exposure
NGROUPS  = int(input('NGROUP')) / number of groups within integration
ZEROFRAM = input('ZROFRAME') / T if a zero frame was read separately
READPATT = input('READOUT') / readout pattern
DATAPROB = F                / T if science telemetry indicated any problems with the data

