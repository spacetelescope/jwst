# Routines used for building cubes

import sys
import time
import numpy as np
import math
import json
import os
from .. import datamodels


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


##################################################################################
class FileTable():
    """
    Dictionary that maps the input files to correct:
    MIRI - channel/subchannel
    NIRSPEC - grating/filter

    Parameters
    ----------

    """
    def __init__(self):

        self.FileMap = {}
        self.FileMap['MIRI'] = {}

        self.FileMap['MIRI']['1'] = {}
        self.FileMap['MIRI']['1']['SHORT'] = []
        self.FileMap['MIRI']['1']['MEDIUM'] = []
        self.FileMap['MIRI']['1']['LONG'] = []

        self.FileMap['MIRI']['2'] = {}
        self.FileMap['MIRI']['2']['SHORT'] = []
        self.FileMap['MIRI']['2']['MEDIUM'] = []
        self.FileMap['MIRI']['2']['LONG'] = []

        self.FileMap['MIRI']['3'] = {}
        self.FileMap['MIRI']['3']['SHORT'] = []
        self.FileMap['MIRI']['3']['MEDIUM'] = []
        self.FileMap['MIRI']['3']['LONG'] = []

        self.FileMap['MIRI']['4'] = {}
        self.FileMap['MIRI']['4']['SHORT'] = []
        self.FileMap['MIRI']['4']['MEDIUM'] = []
        self.FileMap['MIRI']['4']['LONG'] = []

        self.FileMap['NIRSPEC'] = {}
        self.FileMap['NIRSPEC']['PRISM'] = {}
        self.FileMap['NIRSPEC']['PRISM']['CLEAR'] = []

        self.FileMap['NIRSPEC']['G140M'] = {}
        self.FileMap['NIRSPEC']['G140M']['F070LP'] = []
        self.FileMap['NIRSPEC']['G140M']['F100LP'] = []

        self.FileMap['NIRSPEC']['G140H'] = {}
        self.FileMap['NIRSPEC']['G140H']['F070LP'] = []
        self.FileMap['NIRSPEC']['G140H']['F100LP'] = []

        self.FileMap['NIRSPEC']['G235M'] = {}
        self.FileMap['NIRSPEC']['G235M']['F170LP'] = []

        self.FileMap['NIRSPEC']['G235H'] = {}
        self.FileMap['NIRSPEC']['G235H']['F170LP'] = []

        self.FileMap['NIRSPEC']['G395M'] = {}
        self.FileMap['NIRSPEC']['G395M']['F290LP'] = []

        self.FileMap['NIRSPEC']['G395H'] = {}
        self.FileMap['NIRSPEC']['G395H']['F290LP'] = []

        self.FileOffset = {}
        self.FileOffset['1'] = {}
        self.FileOffset['1']['SHORT'] = {}
        self.FileOffset['1']['SHORT']['C1'] = []
        self.FileOffset['1']['SHORT']['C2'] = []
        self.FileOffset['1']['MEDIUM'] = {}
        self.FileOffset['1']['MEDIUM']['C1'] = []
        self.FileOffset['1']['MEDIUM']['C2'] = []
        self.FileOffset['1']['LONG'] = {}
        self.FileOffset['1']['LONG']['C1'] = []
        self.FileOffset['1']['LONG']['C2'] = []

        self.FileOffset['2'] = {}
        self.FileOffset['2']['SHORT'] = {}
        self.FileOffset['2']['SHORT']['C1'] = []
        self.FileOffset['2']['SHORT']['C2'] = []
        self.FileOffset['2']['MEDIUM'] = {}
        self.FileOffset['2']['MEDIUM']['C1'] = []
        self.FileOffset['2']['MEDIUM']['C2'] = []
        self.FileOffset['2']['LONG'] = {}
        self.FileOffset['2']['LONG']['C1'] = []
        self.FileOffset['2']['LONG']['C2'] = []

        self.FileOffset['3'] = {}
        self.FileOffset['3']['SHORT'] = {}
        self.FileOffset['3']['SHORT']['C1'] = []
        self.FileOffset['3']['SHORT']['C2'] = []
        self.FileOffset['3']['MEDIUM'] = {}
        self.FileOffset['3']['MEDIUM']['C1'] = []
        self.FileOffset['3']['MEDIUM']['C2'] = []
        self.FileOffset['3']['LONG'] = {}
        self.FileOffset['3']['LONG']['C1'] = []
        self.FileOffset['3']['LONG']['C2'] = []

        self.FileOffset['4'] = {}
        self.FileOffset['4']['SHORT'] = {}
        self.FileOffset['4']['SHORT']['C1'] = []
        self.FileOffset['4']['SHORT']['C2'] = []
        self.FileOffset['4']['MEDIUM'] = {}
        self.FileOffset['4']['MEDIUM']['C1'] = []
        self.FileOffset['4']['MEDIUM']['C2'] = []
        self.FileOffset['4']['LONG'] = {}
        self.FileOffset['4']['LONG']['C1'] = []
        self.FileOffset['4']['LONG']['C2'] = []

        self.FileOffset['PRISM'] = {}
        self.FileOffset['PRISM']['CLEAR'] = {}
        self.FileOffset['PRISM']['CLEAR']['C1'] = []
        self.FileOffset['PRISM']['CLEAR']['C2'] = []

        self.FileOffset['G140M'] = {}
        self.FileOffset['G140M']['F070LP'] = {}
        self.FileOffset['G140M']['F070LP']['C1'] = []
        self.FileOffset['G140M']['F070LP']['C2'] = []
        self.FileOffset['G140M']['F100LP'] = {}
        self.FileOffset['G140M']['F100LP']['C1'] = []
        self.FileOffset['G140M']['F100LP']['C2'] = []

        self.FileOffset['G140H'] = {}
        self.FileOffset['G140H']['F070LP'] = {}
        self.FileOffset['G140H']['F070LP']['C1'] = []
        self.FileOffset['G140H']['F070LP']['C2'] = []
        self.FileOffset['G140H']['F100LP'] = {}
        self.FileOffset['G140H']['F100LP']['C1'] = []
        self.FileOffset['G140H']['F100LP']['C2'] = []

        self.FileOffset['G235M'] = {}
        self.FileOffset['G235M']['F170LP'] = {}
        self.FileOffset['G235M']['F170LP']['C1'] = []
        self.FileOffset['G235M']['F170LP']['C2'] = []

        self.FileOffset['G235H'] = {}
        self.FileOffset['G235H']['F170LP'] = {}
        self.FileOffset['G235H']['F170LP']['C1'] = []
        self.FileOffset['G235H']['F170LP']['C2'] = []

        self.FileOffset['G395M'] = {}
        self.FileOffset['G395M']['F290LP'] = {}
        self.FileOffset['G395M']['F290LP']['C1'] = []
        self.FileOffset['G395M']['F290LP']['C2'] = []

        self.FileOffset['G395H'] = {}
        self.FileOffset['G395H']['F290LP'] = {}
        self.FileOffset['G395H']['F290LP']['C1'] = []
        self.FileOffset['G395H']['F290LP']['C2'] = []


#********************************************************************************
    def set_file_table(self,
                       input_models,
                       input_filenames,
                       input_ra_offset,
                       input_dec_offset):
#********************************************************************************
        """
        Short Summary
        -------------
        Fill in the MasterTable which holds the files that the cube will be constructed
        from. Since MIRI has 2 channels per image this MASTERTable helps to figure out
        which data needs to be use.
        THe MasterTable for MIRI is broken down by channel and subchannel.
        For each channel/subchannel combination - a file is listed that covers those options
        For NIRSPEC the table contains the Grating and Filter for each file.

        If there is a dither offet file then the master table also holds the ra,dec offset for
        each file.


        Returns
        -------
        MasterTable filled in with files needed
        num: number of files to create cube from
        detector

        """
        num = 0
        num = len(input_filenames)

#________________________________________________________________________________
# Loop over input list of files and assign fill in the MasterTable with filename
# for the correct (channel-subchannel) or (grating-subchannel)
        for i in range(num):

            ifile = input_filenames[i]
            input = input_models[i]

        # Open the input data model & Fill in the FileMap information

            with datamodels.IFUImageModel(input) as input_model:

                detector = input_model.meta.instrument.detector
                instrument = input_model.meta.instrument.name
                assign_wcs = input_model.meta.cal_step.assign_wcs

                if(assign_wcs != 'COMPLETE'):
                    raise ErrorNoAssignWCS("Assign WCS has not been run on file %s",
                                           ifile)
            #________________________________________________________________________________
            #MIRI instrument
            #________________________________________________________________________________
                if instrument == 'MIRI':
                    channel = input_model.meta.instrument.channel
                    subchannel = input_model.meta.instrument.band
            #________________________________________________________________________________
                    clenf = len(channel)
                    for k in range(clenf):
                        self.FileMap['MIRI'][channel[k]][subchannel].append(input_model)
                        ioffset = len(input_ra_offset)
                        if (ioffset > 0):
                            ra_offset = input_ra_offset[i]
                            dec_offset = input_dec_offset[i]
                            self.FileOffset[channel[k]][subchannel]['C1'].append(ra_offset)
                            self.FileOffset[channel[k]][subchannel]['C2'].append(dec_offset)
            #________________________________________________________________________________
                elif instrument== 'NIRSPEC':
                    fwa = input_model.meta.instrument.filter
                    gwa = input_model.meta.instrument.grating

                    self.FileMap['NIRSPEC'][gwa][fwa].append(input_model)
                else:

                    log.info('Instrument not valid for cube')

        return instrument,detector

class ErrorNoAssignWCS(Exception):
    pass
