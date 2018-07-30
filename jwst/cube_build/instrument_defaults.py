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


class InstrumentInfo():

    def __init__(self):
# Wavelength varying parameters

        self.multich_wavelength= []
        self.multich_sroi = []
        self.multich_wroi = []
        self.multich_power = []
        self.multich_softrad = []

        self.prism_wavelength= []
        self.prism_sroi = []
        self.prism_wroi = []
        self.prism_power = []
        self.prism_softrad = []

        self.med_wavelength= []
        self.med_sroi = []
        self.med_wroi = []
        self.med_power = []
        self.med_softrad = []

        self.high_wavelength= []
        self.high_sroi = []
        self.high_wroi = []
        self.high_power = []
        self.high_softrad = []
        
#_______________________________________________________________________
        # This is basic information on the MIRI channels
        # information that will not change (number of slices, starting slice number, ending # slice number and default scales  )
        self.Info = {}
        self.Info['psf_alpha_cuttoff'] = None
        self.Info['psf_alpha_a_short'] = None
        self.Info['psf_alpha_b_short'] = None
        self.Info['psf_alpha_a_long'] = None
        self.Info['psf_alpha_b_long'] = None
        self.Info['psf_beta_cuttoff'] = None
        self.Info['psf_beta_a_short'] = None
        self.Info['psf_beta_b_short'] = None
        self.Info['psf_beta_a_long'] = None
        self.Info['psf_beta_b_long'] = None


#________________________________________________________________________________
        
#________________________________________________________________________________
# channel 1 parameters
        self.Info['1'] = {}
        self.Info['1']['nslices'] = 21
        self.Info['1']['start_slice'] = 101
        self.Info['1']['end_slice'] = 121
        self.Info['1']['xstart'] = 0
        self.Info['1']['xend'] = 512

        self.Info['1']['SHORT'] = {}
        self.Info['1']['SHORT']['ascale'] = None
        self.Info['1']['SHORT']['bscale'] = None
        self.Info['1']['SHORT']['wscale'] = None
        self.Info['1']['SHORT']['sroi'] = None
        self.Info['1']['SHORT']['wroi'] = None
        self.Info['1']['SHORT']['wavemin'] = None
        self.Info['1']['SHORT']['wavemax'] = None
        self.Info['1']['SHORT']['softrad'] = None
        self.Info['1']['SHORT']['msm_power'] = None
        self.Info['1']['SHORT']['rp_wave_cuttoff'] = None
        self.Info['1']['SHORT']['rp_a_low'] = None
        self.Info['1']['SHORT']['rp_b_low'] = None
        self.Info['1']['SHORT']['rp_c_low'] = None
        self.Info['1']['SHORT']['rp_a_high'] = None
        self.Info['1']['SHORT']['rp_b_high'] = None
        self.Info['1']['SHORT']['rp_c_high'] = None
        self.Info['1']['SHORT']['rp_a_ave'] = None
        self.Info['1']['SHORT']['rp_b_ave'] = None
        self.Info['1']['SHORT']['rp_c_ave'] = None

        self.Info['1']['MEDIUM'] = {}
        self.Info['1']['MEDIUM']['ascale'] = None
        self.Info['1']['MEDIUM']['bscale'] = None
        self.Info['1']['MEDIUM']['wscale'] = None
        self.Info['1']['MEDIUM']['sroi'] = None
        self.Info['1']['MEDIUM']['wroi'] = None
        self.Info['1']['MEDIUM']['wavemin'] = None
        self.Info['1']['MEDIUM']['wavemax'] = None
        self.Info['1']['MEDIUM']['softrad'] = None
        self.Info['1']['MEDIUM']['msm_power'] = None
        self.Info['1']['MEDIUM']['rp_wave_cuttoff'] = None
        self.Info['1']['MEDIUM']['rp_a_low'] = None
        self.Info['1']['MEDIUM']['rp_b_low'] = None
        self.Info['1']['MEDIUM']['rp_c_low'] = None
        self.Info['1']['MEDIUM']['rp_a_high'] = None
        self.Info['1']['MEDIUM']['rp_b_high'] = None
        self.Info['1']['MEDIUM']['rp_c_high'] = None
        self.Info['1']['MEDIUM']['rp_a_ave'] = None
        self.Info['1']['MEDIUM']['rp_b_ave'] = None
        self.Info['1']['MEDIUM']['rp_c_ave'] = None

        self.Info['1']['LONG'] = {}
        self.Info['1']['LONG']['ascale'] = None
        self.Info['1']['LONG']['bscale'] = None
        self.Info['1']['LONG']['wscale'] = None
        self.Info['1']['LONG']['sroi'] = None
        self.Info['1']['LONG']['wroi'] = None
        self.Info['1']['LONG']['wavemin'] = None
        self.Info['1']['LONG']['wavemax'] = None
        self.Info['1']['LONG']['softrad'] = None
        self.Info['1']['LONG']['msm_power'] = None
        self.Info['1']['LONG']['rp_wave_cuttoff'] = None
        self.Info['1']['LONG']['rp_a_low'] = None
        self.Info['1']['LONG']['rp_b_low'] = None
        self.Info['1']['LONG']['rp_c_low'] = None
        self.Info['1']['LONG']['rp_a_high'] = None
        self.Info['1']['LONG']['rp_b_high'] = None
        self.Info['1']['LONG']['rp_c_high'] = None
        self.Info['1']['LONG']['rp_a_ave'] = None
        self.Info['1']['LONG']['rp_b_ave'] = None
        self.Info['1']['LONG']['rp_c_ave'] = None
#________________________________________________________________________________
# channel 2 parameters
        self.Info['2'] = {}
        self.Info['2']['nslices'] = 17
        self.Info['2']['start_slice'] = 201
        self.Info['2']['end_slice'] = 217
        self.Info['2']['xstart'] = 513
        self.Info['2']['xend'] = 1031

        self.Info['2']['SHORT'] = {}
        self.Info['2']['SHORT']['ascale'] = None
        self.Info['2']['SHORT']['bscale'] = None
        self.Info['2']['SHORT']['wscale'] = None
        self.Info['2']['SHORT']['sroi'] = None
        self.Info['2']['SHORT']['wroi'] = None
        self.Info['2']['SHORT']['wavemin'] = None
        self.Info['2']['SHORT']['wavemax'] = None
        self.Info['2']['SHORT']['softrad'] = None
        self.Info['2']['SHORT']['msm_power'] = None
        self.Info['2']['SHORT']['rp_wave_cuttoff'] = None
        self.Info['2']['SHORT']['rp_a_low'] = None
        self.Info['2']['SHORT']['rp_b_low'] = None
        self.Info['2']['SHORT']['rp_c_low'] = None
        self.Info['2']['SHORT']['rp_a_high'] = None
        self.Info['2']['SHORT']['rp_b_high'] = None
        self.Info['2']['SHORT']['rp_c_high'] = None
        self.Info['2']['SHORT']['rp_a_ave'] = None
        self.Info['2']['SHORT']['rp_b_ave'] = None
        self.Info['2']['SHORT']['rp_c_ave'] = None

        self.Info['2']['MEDIUM'] = {}
        self.Info['2']['MEDIUM']['ascale'] = None
        self.Info['2']['MEDIUM']['bscale'] = None
        self.Info['2']['MEDIUM']['wscale'] = None
        self.Info['2']['MEDIUM']['sroi'] = None
        self.Info['2']['MEDIUM']['wroi'] = None
        self.Info['2']['MEDIUM']['wavemin'] = None
        self.Info['2']['MEDIUM']['wavemax'] = None
        self.Info['2']['MEDIUM']['softrad'] = None
        self.Info['2']['MEDIUM']['msm_power'] = None
        self.Info['2']['MEDIUM']['rp_wave_cuttoff'] = None
        self.Info['2']['MEDIUM']['rp_a_low'] = None
        self.Info['2']['MEDIUM']['rp_b_low'] = None
        self.Info['2']['MEDIUM']['rp_c_low'] = None
        self.Info['2']['MEDIUM']['rp_a_high'] = None
        self.Info['2']['MEDIUM']['rp_b_high'] = None
        self.Info['2']['MEDIUM']['rp_c_high'] = None
        self.Info['2']['MEDIUM']['rp_a_ave'] = None
        self.Info['2']['MEDIUM']['rp_b_ave'] = None
        self.Info['2']['MEDIUM']['rp_c_ave'] = None

        self.Info['2']['LONG'] = {}
        self.Info['2']['LONG']['ascale'] = None
        self.Info['2']['LONG']['bscale'] = None
        self.Info['2']['LONG']['wscale'] = None
        self.Info['2']['LONG']['sroi'] = None
        self.Info['2']['LONG']['wroi'] = None
        self.Info['2']['LONG']['wavemin'] = None
        self.Info['2']['LONG']['wavemax'] = None
        self.Info['2']['LONG']['softrad'] = None
        self.Info['2']['LONG']['msm_power'] = None
        self.Info['2']['LONG']['rp_wave_cuttoff'] = None
        self.Info['2']['LONG']['rp_a_low'] = None
        self.Info['2']['LONG']['rp_b_low'] = None
        self.Info['2']['LONG']['rp_c_low'] = None
        self.Info['2']['LONG']['rp_a_high'] = None
        self.Info['2']['LONG']['rp_b_high'] = None
        self.Info['2']['LONG']['rp_c_high'] = None
        self.Info['2']['LONG']['rp_a_ave'] = None
        self.Info['2']['LONG']['rp_b_ave'] = None
        self.Info['2']['LONG']['rp_c_ave'] = None
#________________________________________________________________________________
# channel 3 parameters
        self.Info['3'] = {}
        self.Info['3']['nslices'] = 16
        self.Info['3']['start_slice'] = 301
        self.Info['3']['end_slice'] = 316
        self.Info['3']['xstart'] = 513
        self.Info['3']['xend'] = 1031

        self.Info['3']['SHORT'] = {}
        self.Info['3']['SHORT']['ascale'] = None
        self.Info['3']['SHORT']['bscale'] = None
        self.Info['3']['SHORT']['wscale'] = None
        self.Info['3']['SHORT']['sroi'] = None
        self.Info['3']['SHORT']['wroi'] = None
        self.Info['3']['SHORT']['wavemin'] = None
        self.Info['3']['SHORT']['wavemax'] = None
        self.Info['3']['SHORT']['softrad'] = None
        self.Info['3']['SHORT']['msm_power'] = None
        self.Info['3']['SHORT']['rp_wave_cuttoff'] = None
        self.Info['3']['SHORT']['rp_a_low'] = None
        self.Info['3']['SHORT']['rp_b_low'] = None
        self.Info['3']['SHORT']['rp_c_low'] = None
        self.Info['3']['SHORT']['rp_a_high'] = None
        self.Info['3']['SHORT']['rp_b_high'] = None
        self.Info['3']['SHORT']['rp_c_high'] = None
        self.Info['3']['SHORT']['rp_a_ave'] = None
        self.Info['3']['SHORT']['rp_b_ave'] = None
        self.Info['3']['SHORT']['rp_c_ave'] = None

        self.Info['3']['MEDIUM'] = {}
        self.Info['3']['MEDIUM']['ascale'] = None
        self.Info['3']['MEDIUM']['bscale'] = None
        self.Info['3']['MEDIUM']['wscale'] = None
        self.Info['3']['MEDIUM']['sroi'] = None
        self.Info['3']['MEDIUM']['wroi'] = None
        self.Info['3']['MEDIUM']['wavemin'] = None
        self.Info['3']['MEDIUM']['wavemax'] = None
        self.Info['3']['MEDIUM']['softrad'] = None
        self.Info['3']['MEDIUM']['msm_power'] = None
        self.Info['3']['MEDIUM']['rp_wave_cuttoff'] = None
        self.Info['3']['MEDIUM']['rp_a_low'] = None
        self.Info['3']['MEDIUM']['rp_b_low'] = None
        self.Info['3']['MEDIUM']['rp_c_low'] = None
        self.Info['3']['MEDIUM']['rp_a_high'] = None
        self.Info['3']['MEDIUM']['rp_b_high'] = None
        self.Info['3']['MEDIUM']['rp_c_high'] = None
        self.Info['3']['MEDIUM']['rp_a_ave'] = None
        self.Info['3']['MEDIUM']['rp_b_ave'] = None
        self.Info['3']['MEDIUM']['rp_c_ave'] = None

        self.Info['3']['LONG'] = {}
        self.Info['3']['LONG']['ascale'] = None
        self.Info['3']['LONG']['bscale'] = None
        self.Info['3']['LONG']['wscale'] = None
        self.Info['3']['LONG']['sroi'] = None
        self.Info['3']['LONG']['wroi'] = None
        self.Info['3']['LONG']['wavemin'] = None
        self.Info['3']['LONG']['wavemax'] = None
        self.Info['3']['LONG']['softrad'] = None
        self.Info['3']['LONG']['msm_power'] = None
        self.Info['3']['LONG']['rp_wave_cuttoff'] = None
        self.Info['3']['LONG']['rp_a_low'] = None
        self.Info['3']['LONG']['rp_b_low'] = None
        self.Info['3']['LONG']['rp_c_low'] = None
        self.Info['3']['LONG']['rp_a_high'] = None
        self.Info['3']['LONG']['rp_b_high'] = None
        self.Info['3']['LONG']['rp_c_high'] = None
        self.Info['3']['LONG']['rp_a_ave'] = None
        self.Info['3']['LONG']['rp_b_ave'] = None
        self.Info['3']['LONG']['rp_c_ave'] = None
#________________________________________________________________________________
# channel 4 parameters
        self.Info['4'] = {}
        self.Info['4']['nslices'] = 12
        self.Info['4']['start_slice'] = 401
        self.Info['4']['end_slice'] = 412
        self.Info['4']['xstart'] = 0
        self.Info['4']['xend'] = 512

        self.Info['4']['SHORT'] = {}
        self.Info['4']['SHORT']['ascale'] = None
        self.Info['4']['SHORT']['bscale'] = None
        self.Info['4']['SHORT']['wscale'] = None
        self.Info['4']['SHORT']['sroi'] = None        
        self.Info['4']['SHORT']['wroi'] = None
        self.Info['4']['SHORT']['wavemin'] = None
        self.Info['4']['SHORT']['wavemax'] = None
        self.Info['4']['SHORT']['softrad'] = None
        self.Info['4']['SHORT']['msm_power'] = None
        self.Info['4']['SHORT']['rp_wave_cuttoff'] = None
        self.Info['4']['SHORT']['rp_a_low'] = None
        self.Info['4']['SHORT']['rp_b_low'] = None
        self.Info['4']['SHORT']['rp_c_low'] = None
        self.Info['4']['SHORT']['rp_a_high'] = None
        self.Info['4']['SHORT']['rp_b_high'] = None
        self.Info['4']['SHORT']['rp_c_high'] = None
        self.Info['4']['SHORT']['rp_a_ave'] = None
        self.Info['4']['SHORT']['rp_b_ave'] = None
        self.Info['4']['SHORT']['rp_c_ave'] = None

        self.Info['4']['MEDIUM'] = {}
        self.Info['4']['MEDIUM']['ascale'] = None
        self.Info['4']['MEDIUM']['bscale'] = None
        self.Info['4']['MEDIUM']['wscale'] = None
        self.Info['4']['MEDIUM']['sroi'] = None
        self.Info['4']['MEDIUM']['wroi'] = None
        self.Info['4']['MEDIUM']['wavemin'] = None
        self.Info['4']['MEDIUM']['wavemax'] = None
        self.Info['4']['MEDIUM']['softrad'] = None
        self.Info['4']['MEDIUM']['msm_power'] = None
        self.Info['4']['MEDIUM']['rp_wave_cuttoff'] = None
        self.Info['4']['MEDIUM']['rp_a_low'] = None
        self.Info['4']['MEDIUM']['rp_b_low'] = None
        self.Info['4']['MEDIUM']['rp_c_low'] = None
        self.Info['4']['MEDIUM']['rp_a_high'] = None
        self.Info['4']['MEDIUM']['rp_b_high'] = None
        self.Info['4']['MEDIUM']['rp_c_high'] = None
        self.Info['4']['MEDIUM']['rp_a_ave'] = None
        self.Info['4']['MEDIUM']['rp_b_ave'] = None
        self.Info['4']['MEDIUM']['rp_c_ave'] = None

        self.Info['4']['LONG'] = {}
        self.Info['4']['LONG']['ascale'] = None
        self.Info['4']['LONG']['bscale'] = None
        self.Info['4']['LONG']['wscale'] = None
        self.Info['4']['LONG']['wroi'] = None
        self.Info['4']['LONG']['sroi'] = None
        self.Info['4']['LONG']['wavemin'] = None
        self.Info['4']['LONG']['wavemax'] = None
        self.Info['4']['LONG']['softrad'] = None
        self.Info['4']['LONG']['msm_power'] = None
        self.Info['4']['LONG']['rp_wave_cuttoff'] = None
        self.Info['4']['LONG']['rp_a_low'] = None
        self.Info['4']['LONG']['rp_b_low'] = None
        self.Info['4']['LONG']['rp_c_low'] = None
        self.Info['4']['LONG']['rp_a_high'] = None
        self.Info['4']['LONG']['rp_b_high'] = None
        self.Info['4']['LONG']['rp_c_high'] = None
        self.Info['4']['LONG']['rp_a_ave'] = None
        self.Info['4']['LONG']['rp_b_ave'] = None
        self.Info['4']['LONG']['rp_c_ave'] = None

#################################################################################
#NIRSPEC Paramters
        self.Info['PRISM'] = {}
        self.Info['PRISM']['CLEAR'] = {}
        self.Info['PRISM']['CLEAR']['nslices'] = 30
        self.Info['PRISM']['CLEAR']['wscale'] = 0.005
        self.Info['PRISM']['CLEAR']['ascale'] = 0.1
        self.Info['PRISM']['CLEAR']['bscale'] = 0.1
        self.Info['PRISM']['CLEAR']['wroi'] = None
        self.Info['PRISM']['CLEAR']['sroi'] = None
        self.Info['PRISM']['CLEAR']['wavemin'] = None
        self.Info['PRISM']['CLEAR']['wavemax'] = None
        self.Info['PRISM']['CLEAR']['softrad'] = None
        self.Info['PRISM']['CLEAR']['msm_power'] = None

        self.Info['G140M'] = {}
        self.Info['G140M']['F070LP'] = {}
        self.Info['G140M']['F070LP']['nslices'] = 30
        self.Info['G140M']['F070LP']['wscale'] = 0.000636
        self.Info['G140M']['F070LP']['ascale'] = 0.1
        self.Info['G140M']['F070LP']['bscale'] = 0.1
        self.Info['G140M']['F070LP']['wroi'] = None
        self.Info['G140M']['F070LP']['sroi'] = None
        self.Info['G140M']['F070LP']['wavemin'] = None
        self.Info['G140M']['F070LP']['wavemax'] = None
        self.Info['G140M']['F070LP']['softrad'] = None
        self.Info['G140M']['F070LP']['msm_power'] = None

        self.Info['G140M']['F100LP'] = {}
        self.Info['G140M']['F100LP']['nslices'] = 30
        self.Info['G140M']['F100LP']['wscale'] = 0.000636
        self.Info['G140M']['F100LP']['ascale'] = 0.1
        self.Info['G140M']['F100LP']['bscale'] = 0.1
        self.Info['G140M']['F100LP']['wroi'] = None
        self.Info['G140M']['F100LP']['sroi'] = None
        self.Info['G140M']['F100LP']['wavemin'] = None
        self.Info['G140M']['F100LP']['wavemax'] = None
        self.Info['G140M']['F100LP']['softrad'] = None
        self.Info['G140M']['F100LP']['msm_power'] = None

        self.Info['G235M'] = {}
        self.Info['G235M']['F170LP'] = {}
        self.Info['G235M']['F170LP']['nslices'] = 30
        self.Info['G235M']['F170LP']['wscale'] = 0.00106
        self.Info['G235M']['F170LP']['ascale'] = 0.1
        self.Info['G235M']['F170LP']['bscale'] = 0.1
        self.Info['G235M']['F170LP']['wroi'] = None
        self.Info['G235M']['F170LP']['sroi'] = None
        self.Info['G235M']['F170LP']['wavemin'] = None
        self.Info['G235M']['F170LP']['wavemax'] = None
        self.Info['G235M']['F170LP']['softrad'] = None
        self.Info['G235M']['F170LP']['msm_power'] = None

        self.Info['G395M'] = {}
        self.Info['G395M']['F290LP'] = {}
        self.Info['G395M']['F290LP']['nslices'] = 30
        self.Info['G395M']['F290LP']['wscale'] = 0.00179
        self.Info['G395M']['F290LP']['ascale'] = 0.1
        self.Info['G395M']['F290LP']['bscale'] = 0.1
        self.Info['G395M']['F290LP']['wroi'] = None
        self.Info['G395M']['F290LP']['sroi'] = None
        self.Info['G395M']['F290LP']['wavemin'] = None
        self.Info['G395M']['F290LP']['wavemax'] = None
        self.Info['G395M']['F290LP']['softrad'] = None
        self.Info['G395M']['F290LP']['msm_power'] = None


        self.Info['G140H'] = {}
        self.Info['G140H']['F070LP'] = {}
        self.Info['G140H']['F070LP']['nslices'] = 30
        self.Info['G140H']['F070LP']['wscale'] = 0.000235
        self.Info['G140H']['F070LP']['ascale'] = 0.1
        self.Info['G140H']['F070LP']['bscale'] = 0.1
        self.Info['G140H']['F070LP']['wroi'] = None
        self.Info['G140H']['F070LP']['sroi'] = None
        self.Info['G140H']['F070LP']['wavemin'] = None
        self.Info['G140H']['F070LP']['wavemax'] = None
        self.Info['G140H']['F070LP']['softrad'] = None
        self.Info['G140H']['F070LP']['msm_power'] = None

        self.Info['G140H']['F100LP'] = {}
        self.Info['G140H']['F100LP']['nslices'] = 30
        self.Info['G140H']['F100LP']['wscale'] = 0.000235
        self.Info['G140H']['F100LP']['ascale'] = 0.1
        self.Info['G140H']['F100LP']['bscale'] = 0.1
        self.Info['G140H']['F100LP']['wroi'] = None
        self.Info['G140H']['F100LP']['sroi'] = None
        self.Info['G140H']['F100LP']['wavemin'] = None
        self.Info['G140H']['F100LP']['wavemax'] = None
        self.Info['G140H']['F100LP']['softrad'] = None
        self.Info['G140H']['F100LP']['msm_power'] = None

        self.Info['G235H'] = {}
        self.Info['G235H']['F170LP'] = {}
        self.Info['G235H']['F170LP']['nslices'] = 30
        self.Info['G235H']['F170LP']['wscale'] = 0.000396
        self.Info['G235H']['F170LP']['ascale'] = 0.1
        self.Info['G235H']['F170LP']['bscale'] = 0.1
        self.Info['G235H']['F170LP']['wroi'] = None
        self.Info['G235H']['F170LP']['sroi'] = None
        self.Info['G235H']['F170LP']['wavemin'] = None
        self.Info['G235H']['F170LP']['wavemax'] = None
        self.Info['G235H']['F170LP']['softrad'] = None
        self.Info['G235H']['F170LP']['msm_power'] = None

        self.Info['G395H'] = {}
        self.Info['G395H']['F290LP'] = {}
        self.Info['G395H']['F290LP']['nslices'] = 30
        self.Info['G395H']['F290LP']['wscale'] = 0.000665
        self.Info['G395H']['F290LP']['ascale'] = 0.1
        self.Info['G395H']['F290LP']['bscale'] = 0.1
        self.Info['G395H']['F290LP']['wroi'] = None
        self.Info['G395H']['F290LP']['sroi'] = None
        self.Info['G395H']['F290LP']['wavemin'] = None
        self.Info['G395H']['F290LP']['wavemax'] = None
        self.Info['G395H']['F290LP']['softrad'] = None
        self.Info['G395H']['F290LP']['msm_power'] = None

#********************************************************************************
# Functions

    def SetMultiChannelTable(self, wave,sroi,wroi,power,softrad):
        
        self.multich_wavelength.append(wave)
        self.multich_sroi.append(sroi)
        self.multich_wroi.append(wroi)
        self.multich_power.append(power)
        self.multich_softrad.append(softrad)

    def SetPrismTable(self, wave,sroi,wroi,power,softrad):
        self.prism_wavelength.append(wave)
        self.prism_sroi.append(sroi)
        self.prism_wroi.append(wroi)
        self.prism_power.append(power)
        self.prism_softrad.append(softrad)

    def SetMedTable(self, wave,sroi,wroi,power,softrad):
        self.med_wavelength.append(wave)
        self.med_sroi.append(sroi)
        self.med_wroi.append(wroi)
        self.med_power.append(power)
        self.med_softrad.append(softrad)

    def SetHighTable(self, wave,sroi,wroi,power,softrad):
        self.high_wavelength.append(wave)
        self.high_sroi.append(sroi)
        self.high_wroi.append(wroi)
        self.high_power.append(power)
        self.high_softrad.append(softrad)

    def SetSpatialSize(self, value,parameter1,parameter2 = None):
        if parameter2 is None: 
            self.Info[parameter1]['ascale'] = value 
            self.Info[parameter1]['bscale'] = value 
        else: 
            self.Info[parameter1][parameter2]['ascale'] = value 
            self.Info[parameter1][parameter2]['bscale'] = value  


    def SetSpectralStep(self, value,parameter1,parameter2 = None):
        if parameter2 is None: # data is NIRSPEC paramters do not vary with filter
            self.Info[parameter1]['wscale'] = value
        else: 
            self.Info[parameter1][parameter2]['wscale'] = value

    def SetWaveMin(self, value,parameter1,parameter2 = None):
        if parameter2 is None: # data is NIRSPEC paramters do not vary with filter
            self.Info[parameter1]['wavemin'] = value
        else: 
            self.Info[parameter1][parameter2]['wavemin'] = value

    def SetWaveMax(self, value,parameter1,parameter2 = None):
        if parameter2 is None: # data is NIRSPEC paramters do not vary with filter
            self.Info[parameter1]['wavemax'] = value
        else: 
            self.Info[parameter1][parameter2]['wavemax'] = value

    def SetSpatialROI(self, value,parameter1,parameter2 = None):
        if parameter2 is None: # data is NIRSPEC paramters do not vary with filter
            self.Info[parameter1]['sroi'] = value
        else: 
            self.Info[parameter1][parameter2]['sroi'] = value

    def SetWaveROI(self, value,parameter1,parameter2 = None):
        if parameter2 is None: # data is NIRSPEC paramters do not vary with filter
            self.Info[parameter1]['wroi'] = value
        else: 
            self.Info[parameter1][parameter2]['wroi'] = value

    def SetMSMPower(self, value,parameter1,parameter2 = None):
        if parameter2 is None: # data is NIRSPEC paramters do not vary with filter
            self.Info[parameter1]['msm_power'] = value
        else: 
            self.Info[parameter1][parameter2]['msm_power'] = value

    def SetSoftRad(self, value,parameter1,parameter2 = None):
        if parameter2 is None: # data is NIRSPEC paramters do not vary with filter
            self.Info[parameter1]['softrad'] = value
        else: 
            self.Info[parameter1][parameter2]['softrad'] = value

    def Set_RP_Wave_Cutoff(self,table_wave_center,this_channel,this_band):
        self.Info[this_channel][this_band]['rp_wave_cuttoff'] = table_wave_center

    def Set_RP_low(self,a,b,c,this_channel,this_band):
        self.Info[this_channel][this_band]['rp_a_low'] = a
        self.Info[this_channel][this_band]['rp_b_low'] = b
        self.Info[this_channel][this_band]['rp_c_low'] = c

    def Set_RP_high(self,a,b,c,this_channel,this_band):
        self.Info[this_channel][this_band]['rp_a_high'] = a
        self.Info[this_channel][this_band]['rp_b_high'] = b
        self.Info[this_channel][this_band]['rp_c_high'] = c

    def Set_RP_ave(self,a,b,c,this_channel,this_band):
        self.Info[this_channel][this_band]['rp_a_ave'] = a
        self.Info[this_channel][this_band]['rp_b_ave'] = b
        self.Info[this_channel][this_band]['rp_c_ave'] = c


    def Set_psf_alpha_parameters(self,cutoff,a_short,b_short,a_long,b_long):
        self.Info['psf_alpha_cuttoff'] = cutoff
        self.Info['psf_alpha_a_short'] = a_short
        self.Info['psf_alpha_b_short'] = b_short
        self.Info['psf_alpha_a_long'] = a_long
        self.Info['psf_alpha_b_long'] = b_long

    def Set_psf_beta_parameters(self,cutoff,a_short,b_short,a_long,b_long):
        self.Info['psf_beta_cuttoff'] = cutoff
        self.Info['psf_beta_a_short'] = a_short
        self.Info['psf_beta_b_short'] = b_short
        self.Info['psf_beta_a_long'] = a_long
        self.Info['psf_beta_b_long'] = b_long


#______________________________________________________________________
    def Get_RP_ave_Wave(self,this_channel,this_band):
        w = self.Info[this_channel][this_band]['rp_wave_cuttoff']
        a = self.Info[this_channel][this_band]['rp_a_ave']
        b = self.Info[this_channel][this_band]['rp_b_ave']
        c = self.Info[this_channel][this_band]['rp_c_ave']
        weight = (w,a,b,c)
        return weight

    def Get_psf_alpha_parameters(self):
        a1 = self.Info['psf_alpha_cuttoff']
        a2 = self.Info['psf_alpha_a_short']
        a3 = self.Info['psf_alpha_b_short']
        a4 = self.Info['psf_alpha_a_long']
        a5 = self.Info['psf_alpha_b_long']
        a_weight = (a1,a2,a3,a4,a5)
        return a_weight

    def Get_psf_beta_parameters(self):
        b1 = self.Info['psf_beta_cuttoff']
        b2 = self.Info['psf_beta_a_short']
        b3 = self.Info['psf_beta_b_short']
        b4 = self.Info['psf_beta_a_long']
        b5 = self.Info['psf_beta_b_long']
        b_weight = (b1,b2,b3,b4,b5)
        return b_weight

    def GetWaveRoi(self,parameter1,parameter2=None):
        if parameter2 is None:
            roiw = self.Info[parameter1]['wroi']
        else: 
            roiw = self.Info[parameter1][parameter2]['wroi']
        return roiw


    def GetSpatialRoi(self,parameter1,parameter2=None):
        if parameter2 is None:
            rois = self.Info[parameter1]['sroi']
        else: 
            rois = self.Info[parameter1][parameter2]['sroi']
        return rois


    def GetScale(self, parameter1,parameter2 = None):
        if parameter2 is None:
            
            scale = (self.Info[parameter1]['ascale'], 
                     self.Info[parameter1]['bscale'], 
                     self.Info[parameter1]['wscale'])
        else: 
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

    def GetNSlice(self, parameter1):
        numslice = self.Info[parameter1]['nslices']
        return numslice
