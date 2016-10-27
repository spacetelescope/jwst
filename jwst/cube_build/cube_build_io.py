# Routines used for building cubes
from __future__ import absolute_import, print_function

import sys
import time
import numpy as np
import math
import json

from astropy.io import fits

from gwcs.utils import _domain_to_bounds
from ..associations import Association
from .. import datamodels
from ..assign_wcs import nirspec
from . import cube


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


#********************************************************************************
# Read in the User input options for Channel, Subchannel, Filter, Grating

def Read_User_Input(self):
    
    ValidChannel = ['1', '2', '3', '4','ALL']
    ValidSubChannel = ['SHORT', 'MEDIUM', 'LONG','ALL']
    ValidFWA = ['F070LP', 'F100LP', 'F100LP', 'F170LP', 
                    'F170LP', 'F290LP', 'F290LP', 'CLEAR']
    ValidGWA = ['G140M', 'G140H', 'G140M', 'G140H', 'G235M', 'G235H', 
                    'G395M', 'G395H', 'PRISM']
    nchannels = len(ValidChannel)
    nsubchannels = len(ValidSubChannel)

    nfilter = len(ValidFWA)
    ngrating = len(ValidGWA)

#________________________________________________________________________________
    # for MIRI we can set the channel

    if self.channel:  # self.channel is false if it is empty 
        
        channellist = self.channel.split(',')
        user_clen = len(channellist)
#        print('user_clen',user_clen)

#        print('self.channel',self.channel,type(self.channel))
#        print('channellist',channellist)

        for j in range(user_clen):
            ch = channellist[j]
            if(user_clen > 1):
                ch = ch.strip('[')
                ch = ch.strip(']')
                ch = ch.strip(' ')
                ch = ch[1:-1]
            ch = str(ch)

            if ch in ValidChannel:
                self.metadata['channel'].append(ch)
#                print('found channel',ch) 
            else:
                raise ErrorInvalidParameter("Invalid Channel %s",ch)
# remove duplicates if needed
        self.metadata['channel'] = list(set(self.metadata['channel']))


#________________________________________________________________________________
    # for MIRI we can set the subchannel 
    if(self.subchannel): # it is not empty the it has been set
        subchannellist = self.subchannel.split(',')
        user_blen = len(subchannellist)
#        print('number subchannels',user_blen)
        for j in range(user_blen):
            b = subchannellist[j]
            if(user_blen > 1) :
                b = b.strip('[')
                b = b.strip(']')
                b = b.strip(' ')
                b = b[1:-1]
            b  = str(b)
#            print('subchannel',b)
            if b in ValidSubChannel:
                self.metadata['subchannel'].append(b)
#                print('found subchannel',b)
            else:
                raise ErrorInvalidParameter("Invalid Subchannel %s",b)
# remove duplicates if needed
        self.metadata['subchannel'] = list(set(self.metadata['subchannel']))

#________________________________________________________________________________

    # for NIRSPEC we can set the filter
    if(self.filter):
        filterlist = self.filter.split(',')
        user_flen = len(filterlist)
        for j in range(user_flen):
            f = filterlist[j]
            if(user_flen > 1) :
                f = f.strip('[')
                f = f.strip(']')
                f = f.strip(' ')
                f = f[1:-1]
            f  = str(f)
#            print('filter',f)

            if f in ValidFWA:
                self.metadata['filter'].append(f)
#                print('found filter',f)
            else:
                raise ErrorInvalidParameter("Invalid Filter %s",f)
# remove duplicates if needed
        self.metadata['filter'] = list(set(self.metadata['filter']))

#________________________________________________________________________________
    # for NIRSPEC we can set the grating 
    if(self.grating):
        gratinglist = self.grating.split(',')
        user_glen = len(gratinglist)
        for j in range(user_glen):

            g = gratinglist[j]
            if(user_glen > 1) :
                g = g.strip('[')
                g = g.strip(']')
                g = g.strip(' ')
                g = g[1:-1]
            g  = str(g)
#            print('grating',g)
            if g in ValidGWA:
                self.metadata['grating'].append(g)
#                print('found grating',g)
            else:
                raise ErrorInvalidParameter("Invalid Grating %s",g)
# remove duplicates if needed
        self.metadata['grating'] = list(set(self.metadata['grating']))


#********************************************************************************
# Read in dither offset file
# For testing this is useful but possibily this might be useful during flight if the
# images need an additional offset applied to them 
def ReadOffSetFile(self):
#    print('Going to read offset list', self.offset_list)
    f = open(self.offset_list, 'r')
    i = 0
    for line in f:
        #print('input offset', line)
        offset_str = line.split()
        offset = [float(xy) for xy in offset_str]
        ra_off = offset[0]
        dec_off = offset[1]

        self.ra_offset.append(ra_off)
        self.dec_offset.append(dec_off)
#        print('Offset in V2,V3 (in arc seconds) for exposure', ra_off, dec_off, i)
        i = i + 1
    f.close()

#********************************************************************************
def DetermineCubeCoverage(self, MasterTable):
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
    fills in self.metadata['subchannel'] and self.metadata['channel']

    """
#________________________________________________________________________________
# IF INSTRUMENT = MIRI
# loop over the file names

    if(self.metadata['instrument'] == 'MIRI'):
        ValidChannel = ['1', '2', '3', '4']
        ValidSubchannel = ['SHORT', 'MEDIUM', 'LONG']

        nchannels = len(ValidChannel)
        nsubchannels = len(ValidSubchannel)
#________________________________________________________________________________
        # for MIRI we can set the channel and subchannel


        user_clen = len(self.metadata['channel'])
        user_slen = len(self.metadata['subchannel'])

        if(user_clen !=0 and user_slen ==0):
            raise ErrorMissingParameter("Channel specified, but Subchannel was not")

        if(user_clen ==0 and user_slen !=0):
            raise ErrorMissingParameter("Subchannel specified, but Channel was not")

        # parameters not set
        if(user_clen == 0 and  user_slen == 0): 
            for i in range(nchannels):
                for j in range(nsubchannels):
                    nfiles = len(MasterTable.FileMap['MIRI'][ValidChannel[i]][ValidSubchannel[j]])
                    if(nfiles > 0):
                        self.metadata['band_channel'].append(ValidChannel[i])
                        self.metadata['band_subchannel'].append(ValidSubchannel[j])

        # parameters set 
        else:
            for i in range(nchannels):
                for j in range(nsubchannels):
                    nfiles = len(MasterTable.FileMap['MIRI'][ValidChannel[i]][ValidSubchannel[j]])

                    if(nfiles > 0):
                        # now check if these options have been set 
                        if(ValidChannel[i] in self.metadata['channel'] and 
                           ValidSubchannel[i] in self.metadata['subchannel']): 
                            self.metadata['band_channel'].append(ValidChannel[i])
                            self.metadata['band_subchannel'].append(ValidSubchannel[j])


        log.info('The desired cubes covers the MIRI Channels: %s', 
                 self.metadata['band_channel'])
        log.info('The desried cubes covers the MIRI subchannels: %s', 
                 self.metadata['band_subchannel'])


        number_channels = len(self.metadata['band_channel'])
        number_subchannels = len(self.metadata['band_subchannel'])

        if(number_channels == 0):
            raise ErrorNoChannels(
                "The cube  does not cover any channels, change parameter channel")
        if(number_subchannels == 0):
            raise ErrorNoSubchannels(
                "The cube  does not cover any subchannels, change parameter subchannel")
        
        self.metadata['num_bands'] = number_channels # which is = number_subchannels
#______________________________________________________________________
    if(self.metadata['instrument'] == 'NIRSPEC'):

        # 1 to 1 mapping VALIDGWA[i] -> VALIDFWA[i]
        ValidGWA = ['G140M', 'G140H', 'G140M', 'G140H', 'G235M', 'G235H', 
                    'G395M', 'G395H', 'PRISM']
        ValidFWA = ['F070LP', 'F070LP', 'F100LP', 'F100LP', 'F170LP', 
                    'F170LP', 'F290LP', 'F290LP', 'CLEAR']


        nbands = len(ValidFWA)
#________________________________________________________________________________
        # check if input filter or grating has been set
#        grating = self.grating.split()
#        user_glen = len(grating)
#        filter = self.filter.split()
#        user_flen = len(filter) 

        user_glen = len(self.metadata['grating'])
        user_flen = len(self.metadata['filter'])

        if(user_glen ==0 and user_flen !=0):
            raise ErrorMissingParameter("Filter specified, but Grating was not")

        if(user_glen !=0 and user_flen ==0):
            raise ErrorMissingParameter("Grating specified, but Filter was not")
        # Grating and Filter not set - read in from files and create a list of all 
        # the filters and grating contained in the files
        if(user_glen == 0 and  user_flen == 0): 
            for i in range(nbands):
#                print('FWA GWA',ValidFWA[i],ValidGWA[i],nbands)

                nfiles = len(MasterTable.FileMap['NIRSPEC'][ValidGWA[i]][ValidFWA[i]])
                if(nfiles > 0):
                    self.metadata['band_grating'].append(ValidGWA[i])
                    self.metadata['band_filter'].append(ValidFWA[i])

        # Both filter and grating input parameter have been set
        # Find the files that have these parameters set 

        else:
            for i in range(nbands):
                nfiles = len(MasterTable.FileMap['NIRSPEC'][ValidGWA[i]][ValidFWA[i]])
                if(nfiles > 0):
                        # now check if THESE Filter and Grating input parameters were set 
                    if(ValidFWA[i] in self.metadata['filter'] and 
                       ValidGWA[i] in self.metadata['grating']): 
                        self.metadata['band_grating'].append(ValidGWA[i])
                        self.metadata['band_filter'].append(ValidFWA[i])


        log.info('The desired cubes covers the NIRSPEC FWA  %s', 
                  self.metadata['band_filter'])
        log.info('The desried cubes covers the NIRSPEC GWA: %s', 
                  self.metadata['band_grating'])

        number_filters = len(self.metadata['band_filter'])
        number_gratings = len(self.metadata['band_grating'])
        
        self.metadata['num_bands'] = number_gratings # which is = number_filters
        if(number_filters == 0):
            raise ErrorNoFilters("The cube  does not cover any filters")
        if(number_gratings == 0):
            raise ErrorNoGratings("The cube  does not cover any gratings")

#        print('Num gratings and filters',number_gratings,number_filters)

#********************************************************************************
def SetFileTable(self, input_table, MasterTable):
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

    Parameters
    ----------
    input_table: Association Table


    Returns
    -------
    MasterTable filled in with files needed
    num: number of files to create cube from
    detector

    """
    detector = ''
    input_filenames = []
    input_models = []
    i = 0
    num = 0
    iproduct = 0 # only one product found in association table

#________________________________________________________________________________
# find out how many files are in the association table or if it is an single file
# store the input_filenames
    if len(input_table.input_models) > 0:  # this is a single file
        input_models.append(input_table.input_models)
        input_filenames.append(input_table.filename)
        num = 1
    else: # read in assoication table 
        for m in input_table.asn_table['products'][iproduct]['members']:
            input_filenames.append(m['expname'])
            i = i + 1

    num = len(input_filenames)
#    print('number of input filenames',num)
#________________________________________________________________________________
# Loop over input list of files and assign fill in the MasterTable with filename
# for the correct (channel-subchannel) or (grating-subchannel)
    for i in range(num):

        ifile = input_filenames[i]
#        print('openning file',ifile)
        # Open the input data model & Fill in the FileMap information

        with datamodels.ImageModel(ifile) as input_model:

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
                    MasterTable.FileMap['MIRI'][channel[k]][subchannel].append(ifile)
                    ioffset = len(self.ra_offset)
                    if (ioffset > 0):
                        ra_offset = self.ra_offset[i]
                        dec_offset = self.dec_offset[i]
                        MasterTable.FileOffset[channel[k]][subchannel]['C1'].append(ra_offset)
                        MasterTable.FileOffset[channel[k]][subchannel]['C2'].append(dec_offset)
            #________________________________________________________________________________
            elif instrument== 'NIRSPEC':
                fwa = input_model.meta.instrument.filter
                gwa = input_model.meta.instrument.grating

                MasterTable.FileMap['NIRSPEC'][gwa][fwa].append(ifile)
            else:

                print('Instrument not valid for cube')

    return num, instrument,detector

#********************************************************************************
def UpdateOutPutName(self):


    if(self.metadata['instrument'] == 'MIRI'):

        channels = list(set(self.metadata['band_channel']))
        number_channels = len(channels)
        ch_name = '_CH'
        for i in range(number_channels):
            ch_name = ch_name + channels[i]
            if i < number_channels-1:
                ch_name = ch_name + '-'

        
        subchannels = list(set(self.metadata['band_subchannel']))
        number_subchannels = len(subchannels)
        b_name = ''
        for i in range(number_subchannels):
            b_name = b_name + subchannels[i]

        newname = self.output_name_base + ch_name+ '-' + b_name +  '_s3d.fits'

    elif(self.metadata['instrument'] == 'NIRSPEC'):


        fg_name = '_'

        for i in range(self.metadata['num_bands']):
            fg_name = fg_name + self.metadata['band_grating'][i] + '-'+ self.metadata['band_filter'][i]
            if(i < self.metadata['num_bands'] -1):
                fg_name = fg_name + '-'
        newname = self.output_name_base + fg_name+ '_s3d.fits'

    print('Output filename',newname)


    return newname




class ErrorNoAssignWCS(Exception):
    pass

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

        self.input = input # keep a record of original input name for later
        self.input_models = []

        if isinstance(input, datamodels.ImageModel):
            print('this is a single file passed as a Model')
            # It's a single image that's been passed in as a model
            self.interpret_image_model(input)
        elif isinstance(input, str):
            try:
                # The name of an association table
                with open(input, 'r') as input_fh:
                    print('read in association table')
                    self.asn_table = Association.load(input_fh)
            except:
                # The name of a single image file
                print(' this is a single file  read in filename')
                self.filename = input  # temp until figure out model.meta.filename
                self.interpret_image_model(datamodels.ImageModel(input))
        else:
            raise TypeError


    def interpret_image_model(self, model):
        """ Interpret image model as single member association data product.
        """

        # An in-memory ImageModel for a single exposure was provided as input
        self.asn_table = self.template
        self.asn_table['target'] = model.meta.target.catalog_name
        self.asn_table['asn_rule'] = 'singleton'
        self.asn_table['asn_type'] = 'singleton'

        self.asn_table['products'][0]['name'] = self.build_product_name(self.filename)
        self.rootname = self.filename[:self.filename.rfind('_')]
        self.asn_table['products'][0]['members'][0]['expname'] = self.filename
        self.input_models.append(model)

    def build_product_name(self, filename):
        indx = filename.rfind('.fits')
        single_product = filename[:indx]
        return single_product

    def get_inputs(self, product=0):
        members = []
        for p in self.asn_table['products'][product]['members']:
            members.append(p['expname'])
        return members
    def get_outputs(self, product=0):
        return self.asn_table['products'][product]['name']


##################################################################################
class FileTable(object):
    # Dictionary that maps the input files to the 
    # MIRI: Channel & Subchannel
    #NIRSPEC: Grating & Filter
    def __init__(self):

        self.FileMap = {}
        self.FileMap['MIRI'] = {}

        self.FileMap['MIRI']['1'] = {}
        self.FileMap['MIRI']['1']['SHORT'] = list()
        self.FileMap['MIRI']['1']['MEDIUM'] = list()
        self.FileMap['MIRI']['1']['LONG'] = list()

        self.FileMap['MIRI']['2'] = {}
        self.FileMap['MIRI']['2']['SHORT'] = list()
        self.FileMap['MIRI']['2']['MEDIUM'] = list()
        self.FileMap['MIRI']['2']['LONG'] = list()

        self.FileMap['MIRI']['3'] = {}
        self.FileMap['MIRI']['3']['SHORT'] = list()
        self.FileMap['MIRI']['3']['MEDIUM'] = list()
        self.FileMap['MIRI']['3']['LONG'] = list()

        self.FileMap['MIRI']['4'] = {}
        self.FileMap['MIRI']['4']['SHORT'] = list()
        self.FileMap['MIRI']['4']['MEDIUM'] = list()
        self.FileMap['MIRI']['4']['LONG'] = list()

        self.FileMap['NIRSPEC'] = {}
        self.FileMap['NIRSPEC']['PRISM'] = {}
        self.FileMap['NIRSPEC']['PRISM']['CLEAR'] = list()

        self.FileMap['NIRSPEC']['G140M'] = {}
        self.FileMap['NIRSPEC']['G140M']['F070LP'] = list()
        self.FileMap['NIRSPEC']['G140M']['F100LP'] = list()

        self.FileMap['NIRSPEC']['G140H'] = {}
        self.FileMap['NIRSPEC']['G140H']['F070LP'] = list()
        self.FileMap['NIRSPEC']['G140H']['F100LP'] = list()

        self.FileMap['NIRSPEC']['G235M'] = {}
        self.FileMap['NIRSPEC']['G235M']['F170LP'] = list()

        self.FileMap['NIRSPEC']['G235H'] = {}
        self.FileMap['NIRSPEC']['G235H']['F170LP'] = list()

        self.FileMap['NIRSPEC']['G395M'] = {}
        self.FileMap['NIRSPEC']['G395M']['F290LP'] = list()

        self.FileMap['NIRSPEC']['G395H'] = {}
        self.FileMap['NIRSPEC']['G395H']['F290LP'] = list()



        self.FileOffset = {}
        self.FileOffset['1'] = {}
        self.FileOffset['1']['SHORT'] = {}
        self.FileOffset['1']['SHORT']['C1'] = list()
        self.FileOffset['1']['SHORT']['C2'] = list()
        self.FileOffset['1']['MEDIUM'] = {}
        self.FileOffset['1']['MEDIUM']['C1'] = list()
        self.FileOffset['1']['MEDIUM']['C2'] = list()
        self.FileOffset['1']['LONG'] = {}
        self.FileOffset['1']['LONG']['C1'] = list()
        self.FileOffset['1']['LONG']['C2'] = list()

        self.FileOffset['2'] = {}
        self.FileOffset['2']['SHORT'] = {}
        self.FileOffset['2']['SHORT']['C1'] = list()
        self.FileOffset['2']['SHORT']['C2'] = list()

        self.FileOffset['2']['MEDIUM'] = {}
        self.FileOffset['2']['MEDIUM']['C1'] = list()
        self.FileOffset['2']['MEDIUM']['C2'] = list()

        self.FileOffset['2']['LONG'] = {}
        self.FileOffset['2']['LONG']['C1'] = list()
        self.FileOffset['2']['LONG']['C2'] = list()


        self.FileOffset['3'] = {}
        self.FileOffset['3']['SHORT'] = {}
        self.FileOffset['3']['SHORT']['C1'] = list()
        self.FileOffset['3']['SHORT']['C2'] = list()

        self.FileOffset['3']['MEDIUM'] = {}
        self.FileOffset['3']['MEDIUM']['C1'] = list()
        self.FileOffset['3']['MEDIUM']['C2'] = list()

        self.FileOffset['3']['LONG'] = {}
        self.FileOffset['3']['LONG']['C1'] = list()
        self.FileOffset['3']['LONG']['C2'] = list()


        self.FileOffset['4'] = {}
        self.FileOffset['4']['SHORT'] = {}
        self.FileOffset['4']['SHORT']['C1'] = list()
        self.FileOffset['4']['SHORT']['C2'] = list()

        self.FileOffset['4']['MEDIUM'] = {}
        self.FileOffset['4']['MEDIUM']['C1'] = list()
        self.FileOffset['4']['MEDIUM']['C2'] = list()

        self.FileOffset['4']['LONG'] = {}
        self.FileOffset['4']['LONG']['C1'] = list()
        self.FileOffset['4']['LONG']['C2'] = list()

        self.FileOffset['PRISM'] = {}
        self.FileOffset['PRISM']['CLEAR'] = {}
        self.FileOffset['PRISM']['CLEAR']['C1'] = list()
        self.FileOffset['PRISM']['CLEAR']['C2'] = list()

        self.FileOffset['G140M'] = {}
        self.FileOffset['G140M']['F070LP'] = {}
        self.FileOffset['G140M']['F070LP']['C1'] = list()
        self.FileOffset['G140M']['F070LP']['C2'] = list()

        self.FileOffset['G140M']['F100LP'] = {}
        self.FileOffset['G140M']['F100LP']['C1'] = list()
        self.FileOffset['G140M']['F100LP']['C2'] = list()

        self.FileOffset['G140H'] = {}
        self.FileOffset['G140H']['F070LP'] = {}
        self.FileOffset['G140H']['F070LP']['C1'] = list()
        self.FileOffset['G140H']['F070LP']['C2'] = list()

        self.FileOffset['G140H']['F100LP'] = {}
        self.FileOffset['G140H']['F100LP']['C1'] = list()
        self.FileOffset['G140H']['F100LP']['C2'] = list()


        self.FileOffset['G235M'] = {}
        self.FileOffset['G235M']['F170LP'] = {}
        self.FileOffset['G235M']['F170LP']['C1'] = list()
        self.FileOffset['G235M']['F170LP']['C2'] = list()

        self.FileOffset['G235H'] = {}
        self.FileOffset['G235H']['F170LP'] = {}
        self.FileOffset['G235H']['F170LP']['C1'] = list()
        self.FileOffset['G235H']['F170LP']['C2'] = list()


        self.FileOffset['G395M'] = {}
        self.FileOffset['G395M']['F290LP'] = {}
        self.FileOffset['G395M']['F290LP']['C1'] = list()
        self.FileOffset['G395M']['F290LP']['C2'] = list()

        self.FileOffset['G395H'] = {}
        self.FileOffset['G395H']['F290LP'] = {}
        self.FileOffset['G395H']['F290LP']['C1'] = list()
        self.FileOffset['G395H']['F290LP']['C2'] = list()

