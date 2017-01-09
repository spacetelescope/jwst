# Instrument class
# Information on MIRI and NIRSPEC.
# Basic information that will not change
# Default sampling to use based on MIRI:Channel,subchannel, NIRSPEC: FWA,GWS

import sys
import numpy as np
import math
import logging
#from jwst import datamodels


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class Info(object):

    def __init__(self):
#_______________________________________________________________________
        # This is basic information on the MIRI channels
        # information that will not change (number of slices, starting slice number, ending # slice number and default scales  )
        self.Info = {}
        self.Info['1'] = {}
        self.Info['1']['nslices'] = 21
        self.Info['1']['start_slice'] = 101
        self.Info['1']['end_slice'] = 121
        self.Info['1']['xstart'] = 0
        self.Info['1']['xend'] = 512

        self.Info['1']['SHORT'] = {}
        self.Info['1']['SHORT']['ascale'] = 0.13
        self.Info['1']['SHORT']['bscale'] = 0.13
        self.Info['1']['SHORT']['wscale'] = 0.00086975

        self.Info['1']['MEDIUM'] = {}
        self.Info['1']['MEDIUM']['ascale'] = 0.17
        self.Info['1']['MEDIUM']['bscale'] = 0.17
        self.Info['1']['MEDIUM']['wscale'] = 0.00086975

        self.Info['1']['LONG'] = {}
        self.Info['1']['LONG']['ascale'] = 0.17
        self.Info['1']['LONG']['bscale'] = 0.17
        self.Info['1']['LONG']['wscale'] = 0.00086975

        self.Info['2'] = {}
        self.Info['2']['nslices'] = 17
        self.Info['2']['start_slice'] = 201
        self.Info['2']['end_slice'] = 217
        self.Info['2']['xstart'] = 513
        self.Info['2']['xend'] = 1031

        self.Info['2']['SHORT'] = {}
        self.Info['2']['SHORT']['ascale'] = 0.13
        self.Info['2']['SHORT']['bscale'] = 0.13
        self.Info['2']['SHORT']['wscale'] = 0.00154972

        self.Info['2']['MEDIUM'] = {}
        self.Info['2']['MEDIUM']['ascale'] = 0.13
        self.Info['2']['MEDIUM']['bscale'] = 0.13
        self.Info['2']['MEDIUM']['wscale'] = 0.00154972

        self.Info['2']['LONG'] = {}
        self.Info['2']['LONG']['ascale'] = 0.17
        self.Info['2']['LONG']['bscale'] = 0.17
        self.Info['2']['LONG']['wscale'] = 0.00154972

        self.Info['3'] = {}
        self.Info['3']['nslices'] = 16
        self.Info['3']['start_slice'] = 301
        self.Info['3']['end_slice'] = 316
        self.Info['3']['xstart'] = 513
        self.Info['3']['xend'] = 1031

        self.Info['3']['SHORT'] = {}
        self.Info['3']['SHORT']['ascale'] = 0.20
        self.Info['3']['SHORT']['bscale'] = 0.20
        self.Info['3']['SHORT']['wscale'] = 0.002

        self.Info['3']['MEDIUM'] = {}
        self.Info['3']['MEDIUM']['ascale'] = 0.20
        self.Info['3']['MEDIUM']['bscale'] = 0.20
        self.Info['3']['MEDIUM']['wscale'] = 0.002

        self.Info['3']['LONG'] = {}
        self.Info['3']['LONG']['ascale'] = 0.25
        self.Info['3']['LONG']['bscale'] = 0.25
        self.Info['3']['LONG']['wscale'] = 0.002

        self.Info['4'] = {}
        self.Info['4']['nslices'] = 12
        self.Info['4']['start_slice'] = 401
        self.Info['4']['end_slice'] = 412
        self.Info['4']['xstart'] = 0
        self.Info['4']['xend'] = 512

        self.Info['4']['SHORT'] = {}
        self.Info['4']['SHORT']['ascale'] = 0.35
        self.Info['4']['SHORT']['bscale'] = 0.35
        self.Info['4']['SHORT']['wscale'] = 0.003

        self.Info['4']['MEDIUM'] = {}
        self.Info['4']['MEDIUM']['ascale'] = 0.30
        self.Info['4']['MEDIUM']['bscale'] = 0.30
        self.Info['4']['MEDIUM']['wscale'] = 0.003

        self.Info['4']['LONG'] = {}
        self.Info['4']['LONG']['ascale'] = 0.35
        self.Info['4']['LONG']['bscale'] = 0.35
        self.Info['4']['LONG']['wscale'] = 0.003

        self.Info['PRISM'] = {}
        self.Info['PRISM']['wscale'] = 0.005
        self.Info['PRISM']['ascale'] = 0.1
        self.Info['PRISM']['bscale'] = 0.1
        self.Info['PRISM']['nslices'] = 30

        self.Info['G140M'] = {}
        self.Info['G140M']['wscale'] = 0.000636
        self.Info['G140M']['ascale'] = 0.1
        self.Info['G140M']['bscale'] = 0.1
        self.Info['G140M']['nslices'] = 30

        self.Info['G235M'] = {}
        self.Info['G235M']['wscale'] = 0.00106
        self.Info['G235M']['ascale'] = 0.1
        self.Info['G235M']['bscale'] = 0.1
        self.Info['G235M']['nslices'] = 30

        self.Info['G395M'] = {}
        self.Info['G395M']['wscale'] = 0.00179
        self.Info['G395M']['ascale'] = 0.1
        self.Info['G395M']['bscale'] = 0.1
        self.Info['G395M']['bscale'] = 0.1
        self.Info['G395M']['nslices'] = 30

        self.Info['G140H'] = {}
        self.Info['G140H']['wscale'] = 0.000235
        self.Info['G140H']['ascale'] = 0.1
        self.Info['G140H']['bscale'] = 0.1
        self.Info['G140H']['nslices'] = 30

        self.Info['G235H'] = {}
        self.Info['G235H']['wscale'] = 0.000396
        self.Info['G235H']['ascale'] = 0.1
        self.Info['G235H']['bscale'] = 0.1
        self.Info['G235H']['nslices'] = 30

        self.Info['G395H'] = {}
        self.Info['G395H']['wscale'] = 0.000665
        self.Info['G395H']['ascale'] = 0.1
        self.Info['G395H']['bscale'] = 0.1
        self.Info['G395H']['nslices'] = 30

    def GetScale(self, parameter1,parameter2):
        scale = (self.Info[parameter1][parameter2]['ascale'], 
                 self.Info[parameter1][parameter2]['bscale'], 
                 self.Info[parameter1][parameter2]['wscale'])
        return scale


    def GetMIRISliceEndPts(self, parameter1):
        slice_xstart = self.Info[parameter1]['xstart']
        slice_xend = self.Info[parameter1]['xend']
        return slice_xstart, slice_xend

    def GetStartSlice(self, parameter1):
        sliceno = self.Info[parameter1]['start_slice']
        return sliceno

    def GetEndSlice(self, parameter1):
        sliceno = self.Info[parameter1]['end_slice']
        return sliceno
