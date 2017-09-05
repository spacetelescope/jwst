 #! /usr/bin/env python

import sys
import time
import math
import json
import os
import numpy as np
from ..stpipe import Step, cmdline
from fitsblender import blendheaders
from .. import datamodels
from . import cube_build
from . import data_types


class CubeBuildStep (Step):
    """
    CubeBuildStep: Creates a 3-D spectral cube from a given association, single model,
    single input file, or model container.
    Input parameters allow the spectral cube to be built from a provided
    channel/subchannel (MIRI) or grating/filer  (NIRSPEC)
    """

    spec = """
         channel = option('1','2','3','4','ALL','all',default='ALL')
         band = option('SHORT','MEDIUM','LONG','ALL','short','medium','long','all',default='ALL')
         grating   = option('PRISIM','G140M','G140H','G235M','G235H',G395M','G395H','ALL','all',default='ALL')
         filter   = option('CLEAR','F100LP','F070LP','F170LP','F290LP','ALL','all',default='ALL')
         scale1 = float(default=0.0)
         scale2 = float(default=0.0)
         scalew = float(default=0.0)
         weighting = option('msm','miripsf','area','MSM','MIRIPSF','AREA',default = 'msm')
         coord_system = option('ra-dec','alpha-beta','ALPHA-BETA',default='ra-dec')
         rois = float(default=0.0)
         roiw = float(default=0.0)
         weight_power = float(default=2.0)
         offset_list = string(default='NA')
         wavemin = float(default=None)
         wavemax = float(default=None)
         xdebug = integer(default=None)
         ydebug = integer(default=None)
         zdebug = integer(default=None)
         single = boolean(default=false)
         output_type = option('band','channel','grating','multi',default='band')
       """
    reference_file_types = ['cubepar','resol']

    def process(self, input):
        self.log.info('Starting IFU Cube Building Step')

#________________________________________________________________________________
# For all parameters convert to a standard format
# Report read in values to screen
#________________________________________________________________________________
        self.subchannel = self.band
        if(not self.subchannel.isupper()): self.subchannel = self.subchannel.upper()
        if(not self.filter.isupper()): self.filter = self.filter.upper()
        if(not self.grating.isupper()): self.grating = self.grating.upper()
        if(not self.coord_system.islower()): self.coord_system = self.coord_system.lower()
        if(not self.weighting.islower()): self.weighting = self.weighting.lower()

        if(self.scale1 != 0.0): self.log.info('Input Scale of axis 1 %f', self.scale1)
        if(self.scale2 != 0.0): self.log.info('Input Scale of axis 2 %f', self.scale2)
        if(self.scalew != 0.0): self.log.info('Input wavelength scale %f  ', self.scalew)
        if(self.offset_list != 'NA'): self.log.info('Offset Dither list %s', self.offset_list)

        if(self.wavemin !=None): self.log.info('Setting Minimum wavelength of spectral cube to: %f',
                                               self.wavemin)
        if(self.wavemax !=None): self.log.info('Setting Maximum wavelength of spectral cube to: %f',
                                               self.wavemax)

        if(self.rois != 0.0): self.log.info('Input Spatial ROI size %f', self.rois)
        if(self.roiw != 0.0): self.log.info('Input Wave ROI size %f', self.roiw)

        self.debug_pixel = 0
        self.spaxel_debug = None
        if(self.xdebug !=None and self.ydebug !=None and self.zdebug !=None):
            self.debug_pixel = 1
            self.log.info('Writing debug information for spaxel %i %i %i',self.xdebug,self.ydebug,
                          self.zdebug)
            self.log.debug('Writing debug information for spaxel %i %i %i',self.xdebug,self.ydebug,
                           self.zdebug)
            self.xdebug = self.xdebug -1
            self.ydebug = self.ydebug -1
            self.zdebug = self.zdebug -1
            self.spaxel_debug = open('cube_spaxel_info.results','w')
            self.spaxel_debug.write('Writing debug information for spaxel %i %i %i' %
                                    (self.xdebug,self.ydebug,self.zdebug) + '\n')

        # valid coord_system:
        # 1. alpha-beta (only valid for MIRI Single Cubes)
        # 2. ra-dec
        self.interpolation = 'pointcloud' # true for self.weighting  = 'msm' or 'miripsf'

        # if the weighting is area then interpolation is area
        if self.weighting == 'area':
            self.interpolation = 'area'
            self.coord_system = 'alpha-beta'

        if self.coord_system == 'alpha-beta':
            self.weighting = 'area'
            self.interpolation = 'area'

        # if interpolation is point cloud then weighting can be
        # 1. MSM: modified shepard method
        # 2. miripsf - weighting for MIRI based on PSF and LSF
        if self.coord_system == 'ra-dec':
            self.interpolation = 'pointcloud'  # can not be area

        self.log.info('Input interpolation: %s', self.interpolation)
        self.log.info('Coordinate system to use: %s', self.coord_system)
        if self.interpolation =='pointcloud':
            self.log.info('Weighting method for point cloud: %s',self.weighting)
            self.log.info('Power Weighting distance : %f',self.weight_power)

        if self.output_type == 'band':
            self.log.info('Output IFUCubes are Single Bands')
            self.log.info('For MIRI these are single channel, single sub-channel IFU Cubes')
            self.log.info('For NIRSPEC these are single grating, single filter IFU Cubes')

        if self.output_type == 'channel' or self.output_type == 'grating':
            self.log.info('Output IFUCubes are %s', self.output_type)
            self.log.info('For MIRI these are single channel and all subchannels in data')
            self.log.info('   unless the the band option is used to select certain subchannels')
            self.log.info('For NIRSPEC these are single grating and all filters in data')
            self.log.info('   unless the the filter option is used to select certain filters')

        if self.output_type == 'multi':
            self.log.info('Output IFUcube are constructed from all the data unless:')
            self.log.info('   channel, band, grating, or filter options are used')

            
        if self.single :
            self.log.info(' Single = true, creating a set of single exposures mapped' +
                          ' to output IFUCube coordinate system')
            self.output_type = 'single'
#________________________________________________________________________________
    # read input parameters - Channel, Band (Subchannel), Grating, Filter
#________________________________________________________________________________
        self.pars_input = {}
        self.pars_input['channel'] = []     # input parameter or determined from reading in files
        self.pars_input['subchannel'] = []  # inputparameter or determined from reading in files

        self.pars_input['filter'] = []   # input parameter
        self.pars_input['grating'] = []  # input parameter
        read_user_input(self)  # see if options channel, band,grating filter are set
                               # if they are filling par_input with values
#________________________________________________________________________________
#data_types: DataTypes: Read in the input data - 4 formats are allowed:
# 1. filename
# 2. single model
# 3. ASN table
# 4. model containter
# figure out what type of data we have an fill in the
# input_table.input_models - which is used in the rest of IFU Cube Building
# We need to do this in cube_build_step because we need to pass the data_model to
# CRDS to figure out what type of reference files to grab (MIRI or NIRSPEC)
#________________________________________________________________________________
        input_table = data_types.DataTypes(input,self.single)
        self.cube_type = input_table.input_type
        self.input_models = input_table.input_models
        self.input_filenames = input_table.filenames
        self.output_name_base = input_table.output_name
        self.data_type = input_table.data_type
#________________________________________________________________________________
# Read in Cube Parameter Reference file
        # identify what reference file has been associated with these input
        par_filename = self.get_reference_file(self.input_models[0], 'cubepar')
 # Check for a valid reference file
        if par_filename == 'N/A':
            self.log.warning('No default cube parameters reference file found')
            return
#________________________________________________________________________________
# If miripsf weight is set then set up reference file
        # identify what reference file has been associated with these inputs
        resol_filename = None
        if(self.weighting == 'miripsf'):
            resol_filename = self.get_reference_file(self.input_models[0], 'resol')

            if resol_filename == 'N/A':
                self.log.warning('No default spectral resolution reference file found')
                self.log.warning('Run again and turn off miripsf')
                return
#________________________________________________________________________________
# shove the input parameters in to pars to pull out in work horse module -
# cube_build.py

        pars = {
            'channel': self.pars_input['channel'],
            'subchannel': self.pars_input['subchannel'],
            'grating': self.pars_input['grating'],
            'filter': self.pars_input['filter'],
            'scale1': self.scale1,
            'scale2': self.scale2,
            'scalew': self.scalew,
            'interpolation': self.interpolation,
            'weighting': self.weighting,
            'weight_power': self.weight_power,
            'coord_system': self.coord_system,
            'rois': self.rois,
            'roiw': self.roiw,
            'single': self.single,
            'wavemin': self.wavemin,
            'wavemax': self.wavemax,
            'xdebug': self.xdebug,
            'ydebug': self.ydebug,
            'zdebug': self.zdebug,
            'debug_pixel': self.debug_pixel,
            'spaxel_debug':self.spaxel_debug,
            'output_file':self.output_file,
            'output_type':self.output_type,
            'offset_list': self.offset_list}

#________________________________________________________________________________
# create an instance of class CubeData

        cubeinfo = cube_build.CubeData(self.cube_type,
                                       self.input_models,
                                       self.input_filenames,
                                       self.output_name_base,
                                       self.data_type,
                                       par_filename,
                                       resol_filename,
                                       **pars)
#________________________________________________________________________________
# read in all the input files, information from cube_pars, read in input data and
# fill in master_table holding what files are associationed with each ch/sub-ch
# or grating/filter

        self.output_file = cubeinfo.setup()

#        print('in cube_build_step',self.output_file)

#________________________________________________________________________________
# find the min & max final coordinates of cube: map each slice to cube
# add any dither offsets, then find the min & max value in each dimension
# Foot print is returned in ra,dec coordinates

        cubeinfo.setup_wcs()

#________________________________________________________________________________
# build the IFU Cube

# If single = True: map each file to output grid and return single mapped file
#to output grid
# this option is used for background matching and outlier rejection

        if self.single:
            self.output_file = None
            result = cubeinfo.build_ifucube_single()
            self.log.info("Number of IFUCube models returned from building single IFUCubes %i ",len(result))

# Else standard IFU cube building
        else:
           result =  cubeinfo.build_ifucube()
           blendheaders.blendheaders(self.output_file,self.input_filenames)


        if(self.debug_pixel ==1):
            self.spaxel_debug.close()
        return result

#********************************************************************************
# Read in the User input options for Channel, Subchannel, Filter, Grating

def read_user_input(self):
    """
    Short Summary
    -------------
    figure out if any of the input paramters channel,band,filter or grating
    have been set. If they have been  check that they are valid and fill in
    input_pars paramters

    Parameters
    ----------
    none

    Returns
    -------
    self.pars_input['channel']
    self.pars_input['sub_channel']
    self.pars_input['grating']
    self.pars_input['filter']

    """
    ValidChannel = ['1', '2', '3', '4','ALL']
    ValidSubChannel = ['SHORT', 'MEDIUM', 'LONG','ALL']
    ValidFWA = ['F070LP', 'F100LP', 'F100LP', 'F170LP',
                    'F170LP', 'F290LP', 'F290LP', 'CLEAR','ALL']
    ValidGWA = ['G140M', 'G140H', 'G140M', 'G140H', 'G235M', 'G235H',
                    'G395M', 'G395H', 'PRISM','ALL']
    nchannels = len(ValidChannel)
    nsubchannels = len(ValidSubChannel)

    nfilter = len(ValidFWA)
    ngrating = len(ValidGWA)

#________________________________________________________________________________
    # for MIRI we can set the channel
# if set to ALL then let the DetermineCubeCoverage figure out the data we have and set
# self.channel to empty
    if self.channel == 'ALL':
        self.channel = ''

    if self.channel:  # self.channel is false if it is empty

        channellist = self.channel.split(',')
        user_clen = len(channellist)

        for j in range(user_clen):
            ch = channellist[j]
            if(user_clen > 1):
                ch = ch.strip('[')
                ch = ch.strip(']')
                ch = ch.strip(' ')
                ch = ch[1:-1]
            ch = str(ch)

            if ch in ValidChannel:
                self.pars_input['channel'].append(ch)
#                print('found channel',ch)
            else:
                raise ErrorInvalidParameter("Invalid Channel %s",ch)
# remove duplicates if needed
        self.pars_input['channel'] = list(set(self.pars_input['channel']))

#________________________________________________________________________________
    # for MIRI we can set the subchannel
# if set to ALL then let the DetermineCubeCoverage figure out the data we have and set
# self.subchannel = empty

    if self.subchannel == 'ALL':
        self.subchannel = ''

    if self.subchannel : #  not empty it has been set
        subchannellist = self.subchannel.split(',')
        user_blen = len(subchannellist)
        for j in range(user_blen):
            b = subchannellist[j]
            if(user_blen > 1) :
                b = b.strip('[')
                b = b.strip(']')
                b = b.strip(' ')
                b = b[1:-1]
            b  = str(b)
            if b in ValidSubChannel:
                self.pars_input['subchannel'].append(b)
            else:
                raise ErrorInvalidParameter("Invalid Subchannel %s",b)
# remove duplicates if needed
        self.pars_input['subchannel'] = list(set(self.pars_input['subchannel']))
#________________________________________________________________________________
    # for NIRSPEC we can set the filter
# if set to ALL then let the DetermineCubeCoverage figure out the data we have and set
# self.filter = empty
    if self.filter == 'ALL':
        self.filter = ''
    if self.filter:
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
            if f in ValidFWA:
                self.pars_input['filter'].append(f)
            else:
                raise ErrorInvalidParameter("Invalid Filter %s", f)
# remove duplicates if needed
        self.pars_input['filter'] = list(set(self.pars_input['filter']))
#________________________________________________________________________________
    # for NIRSPEC we can set the grating
# if set to ALL then let the DetermineCubeCoverage figure out the data we have and set
# self.grating = empty
    if self.grating == 'ALL':
        self.grating = ''

    if self.grating:
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
            if g in ValidGWA:
                self.pars_input['grating'].append(g)
            else:
                raise ErrorInvalidParameter("Invalid Grating %s",g)
# remove duplicates if needed
        self.pars_input['grating'] = list(set(self.pars_input['grating']))




if __name__ == '__main__':
    cmdline.step_script( cube_build_step )
