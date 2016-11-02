/ Exposure parameters

DETECTOR = input() / name of detector used to acquire data
NINTS    = input() or input('NINT') or 1 / number of integrations within exposure
NGROUPS  = input() or input('NGROUP') / number of groups within integration
ZEROFRAM = input('ZROFRAME') or F  / T if a zero frame was read separately
READPATT = input('READOUT') or 'RAPID' / readout pattern
DATAPROB = F / T if science telemetry indicated any problems

