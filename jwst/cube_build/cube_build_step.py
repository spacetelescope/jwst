 #! /usr/bin/env python

import sys
import time
import math
import json
import os
import numpy as np
from ..stpipe import Step, cmdline
import logging
from fitsblender import blendheaders
from .. import datamodels
from . import cube_build_io
from . import cube_build
from . import cube
from . import CubeCloud
from . import cube_model
from . import InstrumentDefaults
from . import coord
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

1
class CubeBuildStep (Step):
    """
    CubeBuildStep: Creates a 3-D spectral cube from a given association or single input file
    The association will contain which channel/subchannel (MIRI) or filter/grating (NIRSPEC)
    the IFU Cube is going to cover. 
    """

    spec = """
         channel = option('1','2','3','4','ALL','all',default='ALL')
         band = option('SHORT','MEDIUM','LONG','ALL','short','medium','long','all',default='ALL')
         grating   = option('PRISIM','G140M','G140H','G235M','G235H',G395M','G395H','ALL','all',default='ALL')
         filter   = option('CLEAR','F100LP','F070LP','F170LP','F290LP','ALL','all',default='ALL')
         scale1 = float(default=0.0)
         scale2 = float(default=0.0)
         scalew = float(default=0.0)
         interpolation = option(,'pointcloud','area','POINTCLOUD','AREA',default='pointcloud')
         weighting = option('standard','miripsf','STANDARD','MIRIPSF',default = 'standard')
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
       """
    # reference_file_types = ['cubepars','resol']

    def process(self, input):
        self.log.info('Starting IFU Cube Building Step')

        self.subchannel = self.band

        if(not self.subchannel.isupper()): self.subchannel = self.subchannel.upper()
        if(not self.filter.isupper()): self.filter = self.filter.upper()
        if(not self.grating.isupper()): self.grating = self.grating.upper()
        if(not self.coord_system.islower()): self.coord_system = self.coord_system.lower()
        if(not self.interpolation.islower()): self.interpolation = self.interpolation.lower()
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

        if self.interpolation == 'area':
            self.coord_system = 'alpha-beta'

        if self.coord_system == 'ra-dec':
            self.interpolation = 'pointcloud'  # can not be area

        self.log.info('Input interpolation: %s', self.interpolation)
        self.log.info('Coordinate system to use: %s', self.coord_system)
        if self.interpolation =='pointcloud':
            self.log.info('Weighting method for point cloud: %s',self.weighting)
            self.log.info('Power Weighting distance : %f',self.weight_power) 
#_________________________________________________________________________________________________
# Set up the IFU cube basic parameters that define a cube
        self.metadata = {}
        self.metadata['instrument'] = ''
        self.metadata['detector'] = ''
        self.metadata['num_bands'] = 0 

        self.metadata['channel'] = []     # input parameter or determined from reading in files
        self.metadata['subchannel'] = []  # inputparameter or determined from reading in files 

        self.metadata['band_channel'] = []     # band channel: 1-1 pairing with band_subchannel
        self.metadata['band_subchannel'] = []  # band subchannel: 1-1 pairing with band_channel
 
        self.metadata['filter'] = []   # input parameter
        self.metadata['grating'] = []  # input parameter

        self.metadata['band_filter'] = []   # band filter: 1-1 pairing with band_grating 
        self.metadata['band_grating'] = []  # band grating: 1-1 pairing with band_filter

        self.metadata['output_name'] = ''
        self.metadata['number_files'] = 0

        # read input user parameters - Channel, Band (Subchannel), Grating, Filter
        cube_build_io.Read_User_Input(self)
#_________________________________________________________________________________________________
        #Read in the input data - either in form of ASSOCIATION table or single filename
        # If a single file - then assocation table format is filled in

        input_table = cube_build_io.IFUCubeASN(input)

        if input_table.asn_table['asn_type'] == 'singleton':
            self.offset_list = 'NA'

        self.log.info('Output Base %s ', input_table.asn_table['products'][0]['name'])

        # Check if there is an offset list (this ra,dec dither offset list will probably
        # only be used in testing) 
        self.ra_offset = []  # units arc seconds
        self.dec_offset = [] # units arc seconds
        if self.offset_list != 'NA':
            self.log('Going to read in dither offset list')
            cube_build_io.ReadOffSetFile(self)

        # Read in the input data (association table or single file)
        # Fill in MasterTable   based on Channel/Subchannel  or filter/grating
        # Also if there is an Offset list - fill in MasterTable.FileOffset
        MasterTable = cube_build_io.FileTable()
        
        num, instrument, detector = cube_build_io.SetFileTable(self, input_table, 
                                                               MasterTable)

        self.metadata['number_files'] = num
        self.metadata['detector'] = detector            
        self.metadata['instrument'] = instrument

        # self.models is filled in by SetFileTable - only used now for getting the
        # correct reference file type 

        # Determine which channels/subchannels or filter/grating cubes will be constructed from.
        # returns self.metadata['subchannel'] and self.metadata['channel']
        # or self.metadata['filter'], self.metadata['grating']

        cube_build_io.DetermineCubeCoverage(self, MasterTable)
        # check on interpolation = area and coord_system=alpha-beta types of cubes 
        # if interpolation = area also checks that the use did not supply a scale2 values (beta dim) 
        cube_build.CheckCubeType(self)

        self.output_name_base = input_table.asn_table['products'][0]['name']
        self.output_name = cube_build_io.UpdateOutPutName(self)
#________________________________________________________________________________

# Cube is an instance of CubeInfo - which holds basic information on Cube
# Set up if we are building a MIRI cube or a NIRSPEC cube

        if instrument == 'MIRI':
            Cube = cube.CubeInfo('MIRI',detector,
                                 self.metadata['band_channel'], 
                                 self.metadata['band_subchannel'], 
                                 self.output_name,
                                 self.coord_system)

        if instrument == 'NIRSPEC':
            Cube = cube.CubeInfo('NIRSPEC',detector,
                                 self.metadata['band_filter'], 
                                 self.metadata['band_grating'], 
                                 self.output_name,
                                 self.coord_system)

#_________________________________________________________________________________________________
        # for now InstrumentDefaults holds defaults on the two instruments
        InstrumentInfo = InstrumentDefaults.Info()
# Read in Cube Parameter Reference file
        # identify what reference file has been associated with these input
        self.par_filename = self.get_reference_file(self.models[0], 'cubepars')
        self.log.info('Reading default cube parameter file %s', self.par_filename)
        
 # Check for a valid reference file
        if self.par_filename == 'N/A':
            self.log.warning('No default cube parameters reference file found')
            return

        # Load the parameter ref file data model
        # fill in the appropriate fields in InstrumentInfo
        # with the cube parameters

        cube_build_io.ReadCubePars(self,Cube,InstrumentInfo)
#--------------------------------------------------------------------------------
# Read in  MIRI Resolution file if weighting = psf
        # identify what reference file has been associated with these inputs

        if(instrument =='MIRI' and self.weighting == 'miripsf'):
            self.resol_filename = self.get_reference_file(self.models[0], 'resol')
            self.log.info('Reading default MIRI cube resolution file %s', self.resol_filename)
            if self.resol_filename == 'N/A':
                self.log.warning('No default spectral resolution reference file found')
                self.log.warning('Run again and turn off miripsf')

        # Load the miri resolution ref file data model - 
        # fill in the appropriate fields in InstrumentInfo
        # with the cube parameters
            cube_build_io.ReadResolutionFile(self,Cube,InstrumentInfo)
#________________________________________________________________________________
        if self.CubeType == 'File' or self.CubeType == 'ASN' : 
            self.log.info('Building Cube %s ', Cube.output_name)
        # Scale is 3 dimensions and is determined from values held in  InstrumentInfo.GetScale
        scale = cube_build.DetermineScale(Cube, InstrumentInfo)
        # if the user has set the scale of output cube use those values instead
        a_scale = scale[0]
        if self.scale1 != 0.0:
            a_scale = self.scale1

        b_scale = scale[1]
        if self.scale2 != 0.0:
            b_scale = self.scale2
        wscale = scale[2]
        # temp fix for large cubes - need to change to variable wavelength scale
        if self.scalew == 0 and self.metadata['num_bands'] > 6:   
                wscale  = wscale*2            
        if self.scalew == 0 and self.metadata['num_bands'] > 9:   
                wscale  = wscale*2            

        if self.scalew != 0.0:
            wscale = self.scalew

        Cube.SetScale(a_scale, b_scale, wscale)
#________________________________________________________________________________
# get the ROI sizes 
        roi = cube_build.DetermineROISize(Cube, InstrumentInfo)
        # if the user has not set the size of the ROI then use defaults in reference
        # parameter file

        if self.roiw == 0.0: self.roiw = roi[0]
        if self.rois == 0.0: self.rois = roi[1]
        if self.interpolation == 'pointcloud':
            self.log.info('Region of interest  %f %f', 
                               self.rois,self.roiw)
#________________________________________________________________________________
# find the min & max final coordinates of cube: map each slice to cube
# add any dither offsets, then find the min & max value in each dimension
# Foot print is returned in ra,dec coordinates


        # if interpolation = area (alpha-beta), then 1 to 1 mapping in beta 
        CubeFootPrint = cube_build.DetermineCubeSize(self, Cube, 
                                                         MasterTable, 
                                                         InstrumentInfo)
#TODO# do we need the 3 lines of code below ?
#        self.scale1 = Cube.Cdelt1
#        self.scale2 = Cube.Cdelt2
#        self.scalew = Cube.Cdelt3
#        print('new scale2',self.scale2)
#####
        # Based on Scaling and Min and Max values determine naxis1, naxis2, naxis3
        # set cube CRVALs, CRPIXs and xyz coords (center  x,y,z vector spaxel centers)
        if(self.coord_system == 'ra-dec'): 
            Cube.SetGeometry(CubeFootPrint)
        else: 
            Cube.SetGeometryAB(CubeFootPrint) # local coordinate system 

        Cube.PrintCubeGeometry(instrument)
        IFUCube = cube_model.SetUpIFUCube(self,Cube)
#________________________________________________________________________________
        # now you have the size of cube - create an instance for each spaxel
        # create an empty spaxel list - this will become a list of Spaxel classses


        spaxel = []
        total_num = Cube.naxis1*Cube.naxis2*Cube.naxis3

        if(self.interpolation == 'pointcloud'):
            for t in range(total_num):
                spaxel.append(cube.Spaxel())
        else:
            for t in range(total_num):
                spaxel.append(cube.SpaxelAB())

        t0 = time.time()
        # now need to loop over every file that covers this channel/subchannel (MIRI) or Grating/filter(NIRSPEC)
        #and map the detector pixels to the cube spaxel.
        if(instrument == 'MIRI'):
            parameter1 = Cube.channel
            parameter2 = Cube.subchannel
        elif(instrument == 'NIRSPEC'):
            parameter1 = Cube.grating
            parameter2 = Cube.filter

        number_bands = len(parameter1)
        for i in range(number_bands):
            this_par1 = parameter1[i]
            this_par2 = parameter2[i]            
            
            self.log.debug("Working on Band defined by:%s %s " ,this_par1,this_par2)

            cube_build.MapDetectorToCube(self, 
                                         this_par1, this_par2, 
                                         Cube, spaxel, 
                                         MasterTable, 
                                         InstrumentInfo)

        t1 = time.time()
        self.log.info("Time Map All slices on Detector to Cube = %.1f.s" % (t1 - t0,))

#_______________________________________________________________________
# Mapped all data to cube or Point Cloud
# now determine Cube Spaxel flux


        t0 = time.time()
        cube_build.FindCubeFlux(self, Cube, spaxel)

        t1 = time.time()
        self.log.info("Time find Cube Flux= %.1f.s" % (t1 - t0,))

#_______________________________________________________________________
# shove Flux and iflux in the  final IFU cube

        result = cube_model.UpdateIFUCube(self, Cube,IFUCube, spaxel)

#_______________________________________________________________________
 
# write out the IFU cube
#        print(self.CubeType) 
        if self.CubeType == 'File' or self.CubeType =='ASN' :
#            print('Default output file name',self.output_file)
            root, ext = os.path.splitext(self.output_file)
            default = root.find('cube_build') # the user has not provided a name
            if(default != -1):
                self.output_file = IFUCube.meta.filename
        
#_______________________________________________________________________
# Blend the headers

        nf = len(Cube.file)
        print(Cube.file)

        for n in range(nf):
            print('File: ',Cube.file[n])

        
        print('Saving IFUCube',self.output_file)
        IFUCube.save(IFUCube.meta.filename)
#_______________________________________________________________________
        IFUCube.close()

        blendheaders.blendheaders(self.output_file,Cube.file)
#        blendheaders.blendheaders(self.output_file,inputs=Cube.file,verbose=True)

        if(self.debug_pixel ==1):
            self.spaxel_debug.close()
        return result



if __name__ == '__main__':
    cmdline.step_script( cube_build_step )
