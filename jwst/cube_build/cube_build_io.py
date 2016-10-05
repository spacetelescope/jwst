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
# Read in the User input options for Channel, Band, Filter, Grating

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
        print('user_clen',user_clen)

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
                print('found channel',ch) 
            else:
                log.error(' Invalid Channel  %s', ch)
# remove duplicates if needed
        self.metadata['channel'] = list(set(self.metadata['channel']))


#________________________________________________________________________________
    # for MIRI we can set the band 
    if(self.band): # it is not empty the it has been set
        bandlist = self.band.split(',')
        user_blen = len(bandlist)
        print('number bands',user_blen)
        for j in range(user_blen):
            b = bandlist[j]
            if(user_blen > 1) :
                b = b.strip('[')
                b = b.strip(']')
                b = b.strip(' ')
                b = b[1:-1]
            b  = str(b)
            print('band',b)
            if b in ValidSubChannel:
                self.metadata['subchannel'].append(b)
                print('found band',b)
            else:
                log.error(' Invalid SubChannel  %s', b)
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
            print('filter',f)

            if f in ValidFWA:
                self.metadata['filter'].append(f)
                print('found filter',f)
            else:
                log.error(' Invalid Filter give %s', f)
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
            print('grating',g)
            if g in ValidGWA:
                self.metadata['grating'].append(g)
                print('found grating',g)
            else:
                log.error(' Invalid grating give %s', g)
# remove duplicates if needed
        self.metadata['grating'] = list(set(self.metadata['grating']))


    #sys.exit('STOP') 
#********************************************************************************
# Read in dither offset file
# For testing this is useful but possibily this might be useful during flight if the
# images need an additional offset applied to them 
def ReadOffSetFile(self):
    print('Going to read offset list', self.offset_list)
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
        print('Offset in V2,V3 (in arc seconds) for exposure', ra_off, dec_off, i)
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
# loop over the file names

    if(self.metadata['instrument'] == 'MIRI'):
        ValidChannel = ['1', '2', '3', '4']
        ValidSubchannel = ['SHORT', 'MEDIUM', 'LONG']

        nchannels = len(ValidChannel)
        nsubchannels = len(ValidSubchannel)
    

#________________________________________________________________________________
        # for MIRI we can set the channel
        channellist = self.channel.split()
        user_clen = len(channellist)

        # the user has given the channel information
        if(user_clen != 0 ): 
            for j in range(user_clen):

                if channellist[j] in ValidChannel:
                    self.metadata['channel'].append(channellist[j])
                else:
                    log.error(' Invalid Channel give %s', channellist[j])

            self.metadata['channel'] = list(set(self.metadata['channel'])) # remove duplicates if needed

        for i in range(nchannels):
            for j in range(nsubchannels):
                nfiles = len(MasterTable.FileMap['MIRI'][ValidChannel[i]][ValidSubchannel[j]])
                if(nfiles > 0):
                    self.metadata['subchannel'].append(ValidSubchannel[j])
                    if(self.input_table_type == 'singleton' and user_clen == 0):
#usually filled in from reading assoication table
                        self.metadata['channel'].append(ValidChannel[i]) 

        self.metadata['subchannel'] = list(set(self.metadata['subchannel']))
        log.info('The desired cubes covers the MIRI Channels: %s', 
                 self.metadata['channel'])
        log.info('The desried cubes covers the MIRI subchannels: %s', 
                 self.metadata['subchannel'])


        number_channels = len(self.metadata['channel'])
        number_subchannels = len(self.metadata['subchannel'])

        if(number_channels == 0):
            raise ErrorNoChannels(
                "The cube  does not cover any channels, change parameter channel")
        if(number_subchannels == 0):
            raise ErrorNoSubchannels(
                "The cube  does not cover any subchannels, change parameter subchannel")
#______________________________________________________________________
    if(self.metadata['instrument'] == 'NIRSPEC'):

        ValidFWA = ['F070LP', 'F070LP', 'F100LP', 'F100LP', 'F170LP', 
                    'F170LP', 'F290LP', 'F290LP', 'CLEAR']
        ValidGWA = ['G140M', 'G140H', 'G140M', 'G140H', 'G235M', 'G235H', 
                    'G395M', 'G395H', 'PRISM']
        ntypes = len(ValidFWA)

        for j in range(ntypes):
            nfiles = len(MasterTable.FileMap['NIRSPEC'][ValidFWA[j]][ValidGWA[j]])
            if(nfiles > 0):
                self.metadata['filter'].append(ValidFWA[j])
                self.metadata['grating'].append(ValidGWA[j])

        self.metadata['filter'] = list(set(self.metadata['filter']))
        self.metadata['grating'] = list(set(self.metadata['grating']))
        log.debug('The desired cubes covers the NIRSPEC FWA  %s', 
                  self.metadata['filter'])
        log.debug('The desried cubes covers the NIRSPEC GWA: %s', 
                  self.metadata['grating'])

        number_filters = len(self.metadata['filter'])
        number_gratings = len(self.metadata['grating'])

        if(number_filters == 0):
            raise ErrorNoChannels("The cube  does not cover any filters")
        if(number_gratings == 0):
            raise ErrorNoSubchannels("The cube  does not cover any gratings")


#********************************************************************************
def FilesinCube(self, input_table, MasterTable):
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

    if len(input_table.input_models) > 0:  # this is a single file
        input_models.append(input_table.input_models)
        input_filenames.append(input_table.filename)
        num = 1
    else:

        #channels = input_table.asn_table['products'][iproduct]['ch']
#        channels = ['1']
#        channellist = list(channels)
#        num_ch = len(channellist)
#        ValidChannel = ['1', '2', '3', '4']

#        for j in range(num_ch):
#            if channellist[j] in ValidChannel:
#                self.metadata['channel'].append(channellist[j])
#            else:
#                log.error(' Invalid Channel %s', channellist[j])

        for m in input_table.asn_table['products'][iproduct]['members']:

            input_filenames.append(m['expname'])

            i = i + 1

    num = len(input_filenames)
    print('number of input filenames')
#________________________________________________________________________________
# Loop over input list of files
    for i in range(num):

        ifile = input_filenames[i]

        # Open the input data model
        # Fill in the FileMap information

        with datamodels.ImageModel(ifile) as input_model:

            detector = input_model.meta.instrument.detector
            instrument = input_model.meta.instrument.name
            assign_wcs = input_model.meta.cal_step.assign_wcs

            if(assign_wcs != 'COMPLETE'):
                raise ErrorNoAssignWCS("Assign WCS has not been run on file %s", ifile)

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

                MasterTable.FileMap['NIRSPEC'][fwa][gwa].append(ifile)
            else:

                print('Instrument not valid for cube')

    return num, instrument,detector

#********************************************************************************
def UpdateOutPutName(self):

    ch_name = self.channel.replace(" ", "")
    newname = self.output_name + '-CH' + ch_name + '_IFUSCube.fits'
    return newname


#********************************************************************************


#********************************************************************************
def WriteCube(self, Cube, spaxel):

#********************************************************************************
    """
    Short Summary
    -------------
    Write the IFU cube to fits file

    Parameters
    ----------
    Cube: holds meta data of cube
    spaxel: list of spaxels in cube


    Returns
    -------
    no return = writes file

    """
    #pull out data into array

    data = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))
    idata = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))

    dq_cube = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))
    err_cube = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))

    icube = 0
    for z in range(Cube.naxis3):
        for y in range(Cube.naxis2):
            for x in range(Cube.naxis1):
                data[z, y, x] = spaxel[icube].flux
                idata[z, y, x] = len(spaxel[icube].ipointcloud)

                icube = icube + 1
    name = Cube.output_name
    new_model = datamodels.IFUCubeModel(data=data, dq=dq_cube, err=err_cube, weightmap=idata)

    new_model.meta.wcsinfo.crval1 = Cube.Crval1
    new_model.meta.wcsinfo.crval2 = Cube.Crval2
    new_model.meta.wcsinfo.crval3 = Cube.Crval3
    new_model.meta.wcsinfo.crpix1 = Cube.Crpix1
    new_model.meta.wcsinfo.crpix2 = Cube.Crpix2
    new_model.meta.wcsinfo.crpix3 = Cube.Crpix3

    new_model.meta.wcsinfo.crdelt1 = Cube.Cdelt1
    new_model.meta.wcsinfo.crdelt2 = Cube.Cdelt2
    new_model.meta.wcsinfo.crdelt3 = Cube.Cdelt3

    new_model.meta.wcsinfo.ctype1 = 'RA---TAN'
    new_model.meta.wcsinfo.ctype2 = 'DEC--TAN'
    new_model.meta.wcsinfo.ctype3 = 'WAVE'

    new_model.meta.wcsinfo.cunit1 = 'DEG'
    new_model.meta.wcsinfo.cunit2 = 'DEG'
    new_model.meta.wcsinfo.cunit3 = 'MICRON'

#    new_model.meta.wcsinfo.waverange_start = 
#    new_model.meta.wcsinfo.waverange_end = 
    new_model.meta.flux_extension = 'SCI'
    new_model.meta.error_extension = 'ERR'
    new_model.meta.dq_extension = 'DQ'
    new_model.meta.weightmap = 'WMAP'
    new_model.error_type = 'ERR'


    new_model.save(name)
    new_model.close()


    log.info('Wrote %s', name)
    #return new_model

#********************************************************************************
class ErrorNoData(Exception):
    pass

class ErrorNoAssignWCS(Exception):
    pass

class ErrorNoChannels(Exception):
    pass

class ErrorNoSubchannels(Exception):
    pass

class ErrorNoCubes(Exception):
    pass

class IncorrectInput(Exception):
    pass

class NoCoordSystem(Exception):
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

        print('input trying to open',input)
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
                print('Going to read in filename')
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
    # MIRI: Channel & Band
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
        self.FileMap['NIRSPEC']['CLEAR'] = {}
        self.FileMap['NIRSPEC']['CLEAR']['PRISM'] = list()

        self.FileMap['NIRSPEC']['F070LP'] = {}
        self.FileMap['NIRSPEC']['F070LP']['G140M'] = list()
        self.FileMap['NIRSPEC']['F070LP']['G140H'] = list()

        self.FileMap['NIRSPEC']['F100LP'] = {}
        self.FileMap['NIRSPEC']['F100LP']['G140M'] = list()
        self.FileMap['NIRSPEC']['F100LP']['G140H'] = list()

        self.FileMap['NIRSPEC']['F170LP'] = {}
        self.FileMap['NIRSPEC']['F170LP']['G235M'] = list()
        self.FileMap['NIRSPEC']['F170LP']['G235H'] = list()

        self.FileMap['NIRSPEC']['F290LP'] = {}
        self.FileMap['NIRSPEC']['F290LP']['G395M'] = list()
        self.FileMap['NIRSPEC']['F290LP']['G395H'] = list()

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

        self.FileOffset['CLEAR'] = {}
        self.FileOffset['CLEAR']['PRISM'] = {}
        self.FileOffset['CLEAR']['PRISM']['C1'] = list()
        self.FileOffset['CLEAR']['PRISM']['C2'] = list()

        self.FileOffset['F070LP'] = {}
        self.FileOffset['F070LP']['G140M'] = {}
        self.FileOffset['F070LP']['G140M']['C1'] = list()
        self.FileOffset['F070LP']['G140M']['C2'] = list()

        self.FileOffset['F070LP']['G140H'] = {}
        self.FileOffset['F070LP']['G140H']['C1'] = list()
        self.FileOffset['F070LP']['G140H']['C2'] = list()

        self.FileOffset['F100LP'] = {}
        self.FileOffset['F100LP']['G140M'] = {}
        self.FileOffset['F100LP']['G140M']['C1'] = list()
        self.FileOffset['F100LP']['G140M']['C2'] = list()
        self.FileOffset['F100LP']['G140H'] = {}
        self.FileOffset['F100LP']['G140H']['C1'] = list()
        self.FileOffset['F100LP']['G140H']['C2'] = list()

        self.FileOffset['F170LP'] = {}
        self.FileOffset['F170LP']['G235M'] = {}
        self.FileOffset['F170LP']['G235M']['C1'] = list()
        self.FileOffset['F170LP']['G235M']['C2'] = list()
        self.FileOffset['F170LP']['G235H'] = {}
        self.FileOffset['F170LP']['G235H']['C1'] = list()
        self.FileOffset['F170LP']['G235H']['C2'] = list()

        self.FileOffset['F290LP'] = {}
        self.FileOffset['F290LP']['G395M'] = {}
        self.FileOffset['F290LP']['G395M']['C1'] = list()
        self.FileOffset['F290LP']['G395M']['C2'] = list()
        self.FileOffset['F290LP']['G395H'] = {}
        self.FileOffset['F290LP']['G395H']['C1'] = list()
        self.FileOffset['F290LP']['G395H']['C2'] = list()
