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
# determine_band_coverage
# check_cube_type
# update_output_name
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
def determine_band_coverage(self, master_table):
#********************************************************************************
    """
    Short Summary
    -------------
    Function to determine which files contain channels and subchannels are used
    in the creation of the cubes.
    For MIRI The channels  to be used are set by the association and the
    subchannels are  determined from the data

    Parameter
    ----------
    self containing user set input parameters:
    self.channel, self.subchannel

    Returns
    -------
    fills in self.band_channel, self.band_subchannel
    self.band_grating, self.band_filter

    """
#________________________________________________________________________________
# IF INSTRUMENT = MIRI
# loop over the file names

    if self.instrument == 'MIRI':
        ValidChannel = ['1', '2', '3', '4']
        ValidSubchannel = ['SHORT', 'MEDIUM', 'LONG']

        nchannels = len(ValidChannel)
        nsubchannels = len(ValidSubchannel)
#________________________________________________________________________________
        # for MIRI we can set the channel and subchannel

        user_clen = len(self.channel)
        user_slen = len(self.subchannel)
#________________________________________________________________________________
        for i in range(nchannels):
            for j in range(nsubchannels):
                nfiles = len(master_table.FileMap['MIRI'][ValidChannel[i]][ValidSubchannel[j]])
                if nfiles > 0 :
#________________________________________________________________________________
        # neither parameters not set
                    if user_clen == 0 and  user_slen == 0:
                        self.band_channel.append(ValidChannel[i])
                        self.band_subchannel.append(ValidSubchannel[j])
#________________________________________________________________________________
# channel was set by user but not sub-channel
                    elif user_clen !=0 and user_slen ==0:
                        # now check if this channel was set by user
                        if (ValidChannel[i] in self.channel ):
                            self.band_channel.append(ValidChannel[i])
                            self.band_subchannel.append(ValidSubchannel[j])
#________________________________________________________________________________
# sub-channel was set by user but not channel
                    elif user_clen ==0 and user_slen !=0:
                        if (ValidSubchannel[j] in self.subchannel):
                            self.band_channel.append(ValidChannel[i])
                            self.band_subchannel.append(ValidSubchannel[j])
#________________________________________________________________________________
# both parameters set
                    else:
                        if (ValidChannel[i] in self.channel and
                           ValidSubchannel[j] in self.subchannel):
                            self.band_channel.append(ValidChannel[i])
                            self.band_subchannel.append(ValidSubchannel[j])

        log.info('The desired cubes covers the MIRI Channels: %s',
                 self.band_channel)
        log.info('The desired cubes covers the MIRI subchannels: %s',
                 self.band_subchannel)

        number_channels = len(self.band_channel)
        number_subchannels = len(self.band_subchannel)

        if number_channels == 0:
            raise ErrorNoChannels(
                "The cube  does not cover any channels, change parameter channel")
        if number_subchannels == 0:
            raise ErrorNoSubchannels(
                "The cube  does not cover any subchannels, change parameter subchannel")

        self.num_bands = number_channels # which is = number_subchannels
#______________________________________________________________________
    if self.instrument == 'NIRSPEC':

        # 1 to 1 mapping VALIDGWA[i] -> VALIDFWA[i]
        ValidGWA = ['G140M', 'G140H', 'G140M', 'G140H', 'G235M', 'G235H',
                    'G395M', 'G395H', 'PRISM']
        ValidFWA = ['F070LP', 'F070LP', 'F100LP', 'F100LP', 'F170LP',
                    'F170LP', 'F290LP', 'F290LP', 'CLEAR']

        nbands = len(ValidFWA)
#________________________________________________________________________________
        # check if input filter or grating has been set
        user_glen = len(self.grating)
        user_flen = len(self.filter)

        if user_glen ==0 and user_flen !=0:
            raise ErrorMissingParameter("Filter specified, but Grating was not")

        if user_glen !=0 and user_flen ==0:
            raise ErrorMissingParameter("Grating specified, but Filter was not")
        # Grating and Filter not set - read in from files and create a list of all
        # the filters and grating contained in the files
        if user_glen == 0 and  user_flen == 0:
            for i in range(nbands):

                nfiles = len(master_table.FileMap['NIRSPEC'][ValidGWA[i]][ValidFWA[i]])
                if nfiles > 0:
                    self.band_grating.append(ValidGWA[i])
                    self.band_filter.append(ValidFWA[i])

        # Both filter and grating input parameter have been set
        # Find the files that have these parameters set

        else:
            for i in range(nbands):
                nfiles = len(master_table.FileMap['NIRSPEC'][ValidGWA[i]][ValidFWA[i]])
                if nfiles > 0:
                        # now check if THESE Filter and Grating input parameters were set
                    if (ValidFWA[i] in self.filter and
                       ValidGWA[i] in self.grating):
                        self.band_grating.append(ValidGWA[i])
                        self.band_filter.append(ValidFWA[i])


        number_filters = len(self.band_filter)
        number_gratings = len(self.band_grating)

        self.num_bands = number_gratings # which is = number_filters
        if number_filters == 0:
            raise ErrorNoFilters("The cube  does not cover any filters")
        if number_gratings == 0:
            raise ErrorNoGratings("The cube  does not cover any gratings")


#********************************************************************************
def check_cube_type(self):

    if(self.interpolation == "area"):
        if(self.number_files > 1):
            raise IncorrectInput("For interpolation = area, only one file can" +
                                 " be used to created the cube")

        if(len(self.band_channel) > 1):
            raise IncorrectInput("For interpolation = area, only a single channel" +
                                 " can be used to created the cube. Use --channel=# option")

        if(self.scale2 !=0):
            raise AreaInterpolation("When using interpolation = area, the output" +
                                    " coordinate system is alpha-beta" +
                                    " The beta dimension (naxis2) has a one to one" +
                                    " mapping between slice and " +
                                    " beta coordinate.")
                                    

    if(self.coord_system == "alpha-beta"):
        if(self.number_files > 1):
            raise IncorrectInput("Cubes built in alpha-beta coordinate system" +
                                 " are built from a single file")

#********************************************************************************

class IncorrectInput(Exception):
    pass

class NoCoordSystem(Exception):
    pass

class ErrorNoIFUData(Exception):
    pass

class AreaInterpolation(Exception):
    pass

#********************************************************************************
def update_output_name(self):

    if self.cube_type == 'Model':
        newname = self.output_name
    else:
        if self.instrument == 'MIRI':
            #channels = list(set(self.metadata['band_channel']))
            # set does not preserve order so when forming name numbers out of order

            channels = []
            for ch in self.band_channel:
                if ch not in channels:
                       channels.append(ch)

            number_channels = len(channels)
            ch_name = '_ch'
            for i in range(number_channels):
                ch_name = ch_name + channels[i]
                if i < number_channels-1:
                    ch_name = ch_name + '-'


            subchannels = list(set(self.band_subchannel))
            number_subchannels = len(subchannels)
            b_name = ''
            for i in range(number_subchannels):
                b_name = b_name + subchannels[i] 
                if(i > 1): b_name = b_name + '-'
            b_name  = b_name.lower()
            newname = self.output_name_base + ch_name+ '-' + b_name +  '_s3d.fits'

        elif self.instrument == 'NIRSPEC':
            fg_name = '_'

            for i in range(self.num_bands):
                fg_name = fg_name + self.band_grating[i] + '-'+ self.band_filter[i]
                if(i < self.num_bands -1):
                    fg_name = fg_name + '-'
            fg_name = fg_name.lower()

            newname = self.output_name_base + fg_name+ '_s3d.fits'
#________________________________________________________________________________
# check and see if one is provided by the user
# self.output_file is automatically filled in by step class
        #print('output file',self.output_file)
        if(self.output_file == None):
            self.output_file = newname
        else: 
            root, ext = os.path.splitext(self.output_file)
            default = root.find('cube_build') # the user has not provided a name
            if(default != -1):
                self.output_file = newname
            else:
                newname = self.output_file


    return newname


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
        number_bands = len(self.band_channel)

        # pull out the channels and subcahnnels that cover the data making up the cube
        for i in range(number_bands):
            this_channel = self.band_channel[i]
            compare_channel = 'CH'+this_channel
            this_sub = self.band_subchannel[i]
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
        number_gratings = len(self.band_grating)

        for i in range(number_gratings):
            this_gwa = self.band_grating[i]
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
        this_channel = self.band_channel[i]
        this_sub = self.band_subchannel[i]    
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

            for model in input:
                self.input_models.append(model)
                self.filenames.append(model.meta.filename)
#            print('number of models',len(self.filenames))

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

