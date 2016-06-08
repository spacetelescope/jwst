/ Observation identifiers

DATE-OBS  = input('DATE-OBS')  / UTC date at start of exposure
TIME-OBS  = input('TIME-OBS')  / UTC time at start of exposure
OBS_ID    = input()        / full programmatic observation identifier
VISIT_ID  = input()        / visit identifier
PROGRAM   = obsid() or input('OBS_ID')[0:5]   / program number
OBSERVTN  = obsid() or input('OBS_ID')[5:8]   / observation number
VISIT     = obsid() or input('OBS_ID')[8:11]  / visit number
VISITGRP  = obsid() or input('OBS_ID')[11:13] / visit group identifier
SEQ_ID    = obsid() or input('OBS_ID')[13:14] / parallel sequence identifier
ACT_ID    = obsid() or input('OBS_ID')[14:16] / activity identifier
EXPOSURE  = obsid() or input('EXP_ID') / exposure request number

