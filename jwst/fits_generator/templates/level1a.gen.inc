/ Level 1A keywords
INSTRUME  = input() / Instrument used to acquire data
DETECTOR  = { 0x1e1: 'NIR', \
              0x1e2: 'NIR', \
              0x1e3: 'NIR', \
              0x1e4: 'NIR', \
              0x1e5: 'NIR', \
              0x1e6: 'NIR', \
              0x1e7: 'NIR', \
              0x1e8: 'NIR', \
              0x1e9: 'NIR', \
              0x1ea: 'NIR', \
              0x1eb: 'NIR', \
              0x1ec: 'NIR', \
              0x1ed: 'MIR', \
              0x1ee: 'MIR', \
              0x1ef: 'MIR', \
              0x1f0: 'NIR', \
              0x1f1: 'NIR', \
              0x1f2: 'NIR' \
            }[input('SCA_ID')] / Detector type
SCA_NAME  = { 0x1e1: 'NRCA1', \
              0x1e2: 'NRCA2', \
              0x1e3: 'NRCA3', \
              0x1e4: 'NRCA4', \
              0x1e5: 'NRCALONG', \
              0x1e6: 'NRCB1', \
              0x1e7: 'NRCB2', \
              0x1e8: 'NRCB3', \
              0x1e9: 'NRCB4', \
              0x1ea: 'NRCBLONG', \
              0x1eb: 'NRS1', \
              0x1ec: 'NRS2', \
              0x1ed: 'MIRIMAGE', \
              0x1ee: 'MIRIFU1', \
              0x1ef: 'MIRIFU2', \
              0x1f0: 'TFI', \
              0x1f1: 'GUIDER1', \
              0x1f2: 'GUIDER2' \
              }[input('SCA_ID')] / SCA Name
OBS_ID    = input() / Observation ID for the image
PROGRAM   = input('OBS_ID')[0:5] / Program ID
OBSERVTN  = input('OBS_ID')[5:8] / Observation ID
VISIT     = input('OBS_ID')[8:11] / Visit ID
SEQGROUP  = input('OBS_ID')[11:13] / Sequence Group ID
PARASEQN  = input('OBS_ID')[13:14] / Parallel Sequence ID
ACTIVITY  = input('OBS_ID')[14:16] / Activity ID
EXPOSURE  = input('OBS_ID')[16:21] / Exposure ID
NINT      = input() / Count of integration IDs
NGROUP    = input() / Count of readout group
END_TIME  = date_and_time_to_cds(input('DATE-END'), input('TIME-END'))
NCOLUMN   = input('NAXIS1') / Component dimension
NROW      = input('NAXIS2') / Count of readout group
DATAMODE  = {0: 'SLOW', \
             1: 'FAST', \
             2: 'FASTGRPAVG', \
             3: 'FASTINTAVG', \
             25: 'DEEP8', \
             26: 'DEEP2', \
             27: 'MEDIUM8', \
             28: 'MEDIUM2', \
             29: 'SHALLOW4', \
             30: 'SHALLOW2', \
             31: 'BRIGHT2', \
             32: 'BRIGHT1', \
             33: 'RAPID', \
             50: 'NRSRAPID', \
             51: 'NRS', \
             52: 'NRSSLOW', \
             75: 'ID', \
             76: 'ACQ1', \
             77: 'ACQ2', \
             78: 'TRACK', \
             79: 'FINEGUIDE', \
             85: 'TFIGRAPID', \
             86: 'TFIG', \
             } [input()] / Data mode
GAPCOUNT  = input('GROUPGAP') / Data gap counter
COMPRSSD  = bool(input('COMPRESS')) / Compressed flag
CMPLTCND  = input('CMPLTCON')
