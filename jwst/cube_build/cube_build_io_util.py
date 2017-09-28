# Routines used for building cubes
from __future__ import absolute_import, print_function

import sys
import time
import numpy as np
import math
import json
import os

from astropy.io import fits
from ..associations import load_asn
from .. import datamodels
from ..assign_wcs import nirspec
from . import instrument_defaults

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#********************************************************************************
# HELPER ROUTINES for CubeData class defined in cube_build.py
# these methods relate to I/O type procedures.  
# read_offset_file
# read_cubepars
# read_resolution_file
#********************************************************************************
# Read in dither offset file
# For testing this is useful but possibily this might be useful during flight if the
# images need an additional offset applied to them
def read_offset_file(self):

    f = open(self.offset_list, 'r')
    i = 0
    for line in f:
        offset_str = line.split()
        offset = [float(xy) for xy in offset_str]
        ra_off = offset[0]
        dec_off = offset[1]

        self.ra_offset.append(ra_off)
        self.dec_offset.append(dec_off)
        i = i + 1
    f.close()



#********************************************************************************
def read_cubepars(self, instrument_info):
#********************************************************************************
    """
    Short Summary
    -------------
    Based on the instrument and channel/subchannels (MIRI) or grating/filter(NIRSPEC)
    that covers the full range of the data, read in the appropriate columns in the
    cube parameter reference file and fill in the cooresponding dicitionary in 
    instrument_info

    Parameters
    ----------
    ptab: cube parameter reference table 
    instrument_info holds the defaults scales for each channel/subchannel

    Returns
    -------
    The correct elements of instrument_info are filled in 

    """
    if self.instrument == 'MIRI':
        ptab = datamodels.MiriIFUCubeParsModel(self.par_filename)
        number_bands = len(self.all_channel)

        # pull out the channels and subcahnnels that cover the data making up the cube
        for i in range(number_bands):
            this_channel = self.all_channel[i]
            compare_channel = 'CH'+this_channel
            this_sub = self.all_subchannel[i]
                # find the table entries for this combination
            for tabdata in ptab.ifucubepars_table:
                table_channel = tabdata['channel']
                table_band = tabdata['band']
                table_plt_scale = tabdata['PLTSCALE']
                table_wresol = tabdata['WRESOL']
                table_sroi = tabdata['SROI']
                table_wroi = tabdata['WROI']
                #match on this_channel and this_sub
                if(compare_channel == table_channel and this_sub == table_band):
                    instrument_info.SetSpatialScale(table_plt_scale,this_channel,this_sub)
                    instrument_info.SetWaveRes(table_wresol,this_channel,this_sub)
                    instrument_info.SetWaveROI(table_wroi,this_channel,this_sub)
                    instrument_info.SetSpatialROI(table_sroi,this_channel,this_sub)
        
    elif self.instrument == 'NIRSPEC':
        ptab = datamodels.NirspecIFUCubeParsModel(self.par_filename)
        number_gratings = len(self.all_grating)

        for i in range(number_gratings):
            this_gwa = self.all_grating[i]
            for tabdata in ptab.ifucubepars_table:
                table_grating = tabdata['grating']
                table_filter = tabdata['filter']
                table_plt_scale = tabdata['PLTSCALE']
                table_wresol = tabdata['WRESOL']
                table_sroi = tabdata['SROI']
                table_wroi = tabdata['WROI']
                if(this_gwa == table_grating):
                    instrument_info.SetSpatialScale(table_plt_scale,this_gwa)
                    instrument_info.SetWaveRes(table_wresol,this_gwa)
                    instrument_info.SetWaveROI(table_wroi,this_gwa)
                    instrument_info.SetSpatialROI(table_sroi,this_gwa)

#_______________________________________________________________________

# Read MIRI Resolution reference file
#********************************************************************************
def read_resolution_file(self,instrument_info):

    ptab = datamodels.MiriResolutionModel(self.resol_filename)
    table_alpha_cutoff = ptab.psf_fwhm_alpha_table['A_CUTOFF']
    table_alpha_a_short = ptab.psf_fwhm_alpha_table['A_A_SHORT']
    table_alpha_b_short = ptab.psf_fwhm_alpha_table['A_B_SHORT']
    table_alpha_a_long = ptab.psf_fwhm_alpha_table['A_A_LONG']
    table_alpha_b_long = ptab.psf_fwhm_alpha_table['A_B_LONG']

    table_beta_cutoff = ptab.psf_fwhm_beta_table['B_CUTOFF']
    table_beta_a_short = ptab.psf_fwhm_beta_table['B_A_SHORT']
    table_beta_b_short = ptab.psf_fwhm_beta_table['B_B_SHORT']
    table_beta_a_long = ptab.psf_fwhm_beta_table['B_A_LONG']
    table_beta_b_long = ptab.psf_fwhm_beta_table['B_B_LONG']
    

    instrument_info.Set_psf_alpha_parameters(table_alpha_cutoff,
                                            table_alpha_a_short,
                                            table_alpha_b_short,
                                            table_alpha_a_long,
                                            table_alpha_b_long)

    instrument_info.Set_psf_beta_parameters(table_beta_cutoff,
                                            table_beta_a_short,
                                            table_beta_b_short,
                                            table_beta_a_long,
                                            table_beta_b_long)
    
    number_bands = len(self.channel)

        # pull out the channels and subcahnnels that cover the data making up the cube
    for i in range(number_bands):
        this_channel = self.all_channel[i]
        this_sub = self.all_subchannel[i]    
        compare_band = this_channel+this_sub
        for tabdata in ptab.resolving_power_table:
            table_sub_band = tabdata['SUB_BAND']
            table_wave_center = tabdata['R_CENTRE']
            table_res_a_low = tabdata['R_A_LOW']
            table_res_b_low = tabdata['R_B_LOW']
            table_res_c_low = tabdata['R_C_LOW']
            table_res_a_high = tabdata['R_A_HIGH']
            table_res_b_high = tabdata['R_B_HIGH']
            table_res_c_high = tabdata['R_C_HIGH']
            table_res_a_ave = tabdata['R_A_AVG']
            table_res_b_ave = tabdata['R_B_AVG']
            table_res_c_ave = tabdata['R_C_AVG']
            #match on this_channel and this_sub
            if compare_band == table_sub_band:
                instrument_info.Set_RP_Wave_Cutoff(table_wave_center,
                                                  this_channel,this_sub)
                instrument_info.Set_RP_low(table_res_a_low,
                                          table_res_b_low,
                                          table_res_c_low,
                                          this_channel,this_sub)

                instrument_info.Set_RP_high(table_res_a_high,
                                          table_res_b_high,
                                          table_res_c_high,
                                          this_channel,this_sub)

                instrument_info.Set_RP_ave(table_res_a_ave,
                                          table_res_b_ave,
                                          table_res_c_ave,
                                          this_channel,this_sub)


class ErrorNoChannels(Exception):
    pass

class ErrorNoSubchannels(Exception):
    pass

class ErrorNoFilters(Exception):
    pass

class ErrorNoGrating(Exception):
    pass

class ErrorInvalidParameter(Exception):
    pass

class ErrorMissingParameter(Exception):
    pass


#********************************************************************************
class IFUCubeASN(object):
#********************************************************************************

    """
    Class to handle reading the input to the processing, which
    can be a single science exposure or an IFU cube association table.
    The input and output member info is loaded into an ASN table model.
    """

    template = {"asn_rule": "",
              "target": "",
              "asn_pool": "",
              "asn_type": "",
              "products": [
                  {"name": "",
                   "members": [
                      {"exptype": "",
                       "expname": ""}
                      ]
                  }
                ]
              }

    def __init__(self, input):

        self.input_models = []
        self.filenames = []
        self.output_name = None
        self.data_type = None # singleton, multi
        self.input_type = None # Model, File, ASN, Container

        # IF a single model or a single file  is passed in then
        # self.filename & self.input_model hold the values for this singe dataset
        self.InputType  = ''
        if isinstance(input, datamodels.ImageModel):
#            print('this is a single file passed as a Model')
            # It's a single image that's been passed in as a model
            # input is a model
            self.filenames.append(input.meta.filename)
            self.input_models.append(input)
            self.input_type = 'Model'
            self.data_type = 'singleton'
            self.output_name = self.build_product_name(self.filenames[0])
        elif isinstance(input,datamodels.ModelContainer):
#            print('this is a model container type')
            self.input_type='Container'
            self.data_type = 'multi'
            with datamodels.ModelContainer(input) as input_model:
                self.output_name =input_model.meta.asn_table.products[0].name

            for i in range(len(input)):
                model = input[i]
                self.input_models.append(model)
                self.filenames.append(model.meta.filename)

        elif isinstance(input, str):
            try:
                # The name of an association table
                # for associations - use Association.load
                # in cube_build_io.SetFileTable - set up:
                # input_model & filename lists
                iproduct = 0 # only one product found in association table
                with open(input, 'r') as input_fh:
#                    print('read in association table')
                    asn_table = load_asn(input_fh)
                    self.input_type = 'ASN'
                    self.data_type = 'multi'
                    self.output_name =  asn_table['products'][0]['name']
                    for m in asn_table['products'][iproduct]['members']:
                        self.filenames.append(m['expname'])
                        self.input_models.append(datamodels.ImageModel(m['expname']))
            except:
                # The name of a single image file
#                print(' this is a single file  read in filename')
                self.input_type = 'File'
                self.data_type = 'singleton'
                self.filenames.append(input)
                self.input_models.append(datamodels.ImageModel(input))
                self.output_name = self.build_product_name(self.filenames[0])

        else:
            raise TypeError


    def build_product_name(self, filename):
        indx = filename.rfind('.fits')
        single_product = filename[:indx]
        return single_product

# TODO:  Routines not used below - saved just in case we need them later - if not
# remove. 

    def interpret_image_model(self, model):
        """ Interpret image model as single member association data product.
            Currently this routien is not used by cube_build - it was left
            if needed down the road
        """

        # An in-memory ImageModel for a single exposure was provided as input
        self.asn_table = self.template
        self.asn_table['target'] = model.meta.target.catalog_name
        self.asn_table['asn_rule'] = 'singleton'
        self.asn_table['asn_type'] = 'singleton'
        self.asn_table['products'][0]['name'] = self.build_product_name(self.filenames[0])
        self.rootname = self.filename[:self.filename.rfind('_')]
        self.asn_table['products'][0]['members'][0]['expname'] = self.filenames[0]

    def get_inputs(self, product=0):
        members = []
        for p in self.asn_table['products'][product]['members']:
            members.append(p['expname'])
        return members
    def get_outputs(self, product=0):
        return self.asn_table['products'][product]['name']

