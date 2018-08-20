 #! /usr/bin/env python
from ..stpipe import Step
from .. import datamodels
from . import cube_build
from . import ifu_cube
from . import data_types
from ..assign_wcs.util import update_s_region_keyword

__all__ = ["CubeBuildStep", "read_user_input"]

class CubeBuildStep (Step):
    """
    CubeBuildStep: Creates a 3-D spectral cube from a given association, single model,
    single input file, or model container.
    Input parameters allow the spectral cube to be built from a provided
    channel/subchannel (MIRI) or grating/filer  (NIRSPEC)
    """

    spec = """
         channel = option('1','2','3','4','ALL','all',default='ALL') # Options: 1,2,3,4,or All
         band = option('SHORT','MEDIUM','LONG','ALL','short','medium','long','all',default='ALL') # Options: \
SHORT,MEDIUM,LONG, or ALL
         grating   = option('PRISM','G140M','G140H','G235M','G235H',G395M','G395H','ALL','all',default='ALL')  # Options: PRISM,G140M,G140H,G235M,G235H,G395M,G395H, or ALL
         filter   = option('CLEAR','F100LP','F070LP','F170LP','F290LP','ALL','all',default='ALL') # Options: CLEAR,F100LP,F070LP,F170LP,F290LP, or ALL
         scale1 = float(default=0.0) # cube sample size to use for axis 1, arc seconds
         scale2 = float(default=0.0) # cube sample size to use for axis 2, arc seconds
         scalew = float(default=0.0) # cube sample size to use for axis 3, microns
         weighting = option('msm','miripsf','area','MSM','MIRIPSF','AREA',default = 'msm') # Type of weighting function,
         coord_system = option('ra-dec','alpha-beta','ALPHA-BETA',default='ra-dec') # Output Coordinate system. Options: ra-dec or alpha-beta
         rois = float(default=0.0) # region of interest spatial size, arc seconds
         roiw = float(default=0.0) # region of interest wavelength size, microns
         weight_power = float(default=2.0) # Weighting option to use for Modified Shepard Method
         offset_list = string(default='NA')  # A file for dithered data containing additional ra and dec offsets
         wavemin = float(default=None)  # Minimum wavelength to be used in the IFUCube
         wavemax = float(default=None)  # Maximum wavelength to be used in the IFUCube
         single = boolean(default=false) # Internal pipeline option used by mrs_imatch and outlier detection
         xdebug = integer(default=None) # debug option, x spaxel value to report information on
         ydebug = integer(default=None) # debug option, y spaxel value to report information on
         zdebug = integer(default=None) # debug option, z spaxel value to report  information on
         output_type = option('band','channel','grating','multi',default='band') # Type IFUcube to create. Options=band,channel,grating,multi
         search_output_file = boolean(default=false)
         output_use_model = boolean(default=true) # Use filenames in the output models
       """
    reference_file_types = ['cubepar', 'resol']

    def process(self, input):
        self.log.info('Starting IFU Cube Building Step')
#________________________________________________________________________________
# For all parameters convert to a standard format
# Report read in values to screen
#________________________________________________________________________________
        self.subchannel = self.band
        # print('self suffix',self.suffix)
        self.suffix = 's3d' # override suffix = cube_build

        if(not self.subchannel.isupper()): self.subchannel = self.subchannel.upper()
        if(not self.filter.isupper()): self.filter = self.filter.upper()
        if(not self.grating.isupper()): self.grating = self.grating.upper()
        if(not self.coord_system.islower()): self.coord_system = self.coord_system.lower()
        if(not self.output_type.islower()): self.output_type = self.output_type.lower()
        if(not self.weighting.islower()): self.weighting = self.weighting.lower()

        if(self.scale1 != 0.0): self.log.info('Input Scale of axis 1 %f', self.scale1)
        if(self.scale2 != 0.0): self.log.info('Input Scale of axis 2 %f', self.scale2)
        if(self.scalew != 0.0): self.log.info('Input wavelength scale %f  ', self.scalew)
        if(self.offset_list != 'NA'): self.log.info('Offset Dither list %s', self.offset_list)

        if(self.wavemin is not None): self.log.info('Setting Minimum wavelength of spectral cube to: %f',
                                               self.wavemin)
        if(self.wavemax is not None): self.log.info('Setting Maximum wavelength of spectral cube to: %f',
                                               self.wavemax)

        if(self.rois != 0.0): self.log.info('Input Spatial ROI size %f', self.rois)
        if(self.roiw != 0.0): self.log.info('Input Wave ROI size %f', self.roiw)

        self.debug_pixel = 0
        self.spaxel_debug = None
        if(self.xdebug is not None and self.ydebug is not None and self.zdebug is not None):
            self.debug_pixel = 1
            self.log.info('Writing debug information for spaxel %i %i %i',
                          self.xdebug,
                          self.ydebug,
                          self.zdebug)
            self.log.debug('Writing debug information for spaxel %i %i %i',
                           self.xdebug,
                           self.ydebug,
                           self.zdebug)
            self.xdebug = self.xdebug - 1
            self.ydebug = self.ydebug - 1
            self.zdebug = self.zdebug - 1
            self.spaxel_debug = open('cube_spaxel_info.results', 'w')
            self.spaxel_debug.write('Writing debug information for spaxel %i %i %i' %
                                    (self.xdebug, self.ydebug, self.zdebug) + '\n')

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
        if self.interpolation == 'pointcloud':
            self.log.info('Weighting method for point cloud: %s', self.weighting)
            self.log.info('Power Weighting distance : %f', self.weight_power)

        if self.single:
            self.output_type = 'single'
            self.log.info('Cube Type: Single cubes ')

#________________________________________________________________________________
    # read input parameters - Channel, Band (Subchannel), Grating, Filter
#________________________________________________________________________________
        self.pars_input = {}
        self.pars_input['channel'] = []     # input parameter or determined from reading in files
        self.pars_input['subchannel'] = []  # inputparameter or determined from reading in files

        self.pars_input['filter'] = []   # input parameter
        self.pars_input['grating'] = []  # input parameter
        read_user_input(self)  # see if options channel, band,grating filter are set
                               # is they are then self.output_type = 'user'
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
# if the user has provided the filename - strip out .fits and pull out the base name
# the cube_build software will attached the needed information on channel, sub-channel
# grating or filter.
#________________________________________________________________________________
        input_table = data_types.DataTypes(input, self.single,
                                           self.output_file,
                                           self.output_dir)

        self.input_models = input_table.input_models
        self.input_filenames = input_table.filenames
        self.output_name_base = input_table.output_name

        self.pipeline = 3
        if self.output_type == 'multi' and len(self.input_filenames) == 1:
            self.pipeline = 2
        if(len(self.input_filenames) == 1): self.offset_list = 'NA'
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
# shove the input parameters in to pars to pull out in general cube_build.py

        pars = {
            'channel': self.pars_input['channel'],
            'subchannel': self.pars_input['subchannel'],
            'grating': self.pars_input['grating'],
            'filter': self.pars_input['filter'],
            'weighting': self.weighting,
            'single': self.single,
            'output_type': self.output_type,
            'offset_list': self.offset_list}

# shove the input parameters in to pars_cube to pull out ifu_cube.py
# these parameters are related to the IFUCube
        pars_cube = {
            'scale1': self.scale1,
            'scale2': self.scale2,
            'scalew': self.scalew,
            'interpolation': self.interpolation,
            'weighting': self.weighting,
            'weight_power': self.weight_power,
            'coord_system': self.coord_system,
            'rois': self.rois,
            'roiw': self.roiw,
            'wavemin': self.wavemin,
            'wavemax': self.wavemax,
            'xdebug': self.xdebug,
            'ydebug': self.ydebug,
            'zdebug': self.zdebug,
            'debug_pixel': self.debug_pixel,
            'spaxel_debug': self.spaxel_debug}
#________________________________________________________________________________
# create an instance of class CubeData

        cubeinfo = cube_build.CubeData(
            self.input_models,
            self.input_filenames,
            par_filename,
            resol_filename,
            **pars)
#________________________________________________________________________________
# cubeinfo.setup:
# read in all the input files, information from cube_pars, read in input data and
# fill in master_table holding what files are associationed with each ch/sub-ch
# or grating/filter -> fills in all_channel, all_subchannel,all_filter, all_grating
# instrument and detector

        result = cubeinfo.setup()
        instrument = result['instrument']
        detector = result['detector']
        instrument_info = result['instrument_info']
        master_table = result['master_table']
#________________________________________________________________________________
# How many and what type of cubes will be made
# send self.output_type, all_channel, all_subchannel, all_grating, all_filter
# return number of cubes and for each cube the list_pars1, list_pars2 (channel,subchannel)
# or (grating,filter)

        num_cubes, cube_pars = cubeinfo.number_cubes()
        self.log.info('Number of IFUCubes produced by a this run %i', num_cubes)

        cube_container = datamodels.ModelContainer() # ModelContainer of ifucubes

        for i in range(num_cubes):
            icube = str(i + 1)
            list_par1 = cube_pars[icube]['par1']
            list_par2 = cube_pars[icube]['par2']
            thiscube = ifu_cube.IFUCubeData(
                self.pipeline,
                self.input_filenames,
                self.input_models,
                self.output_name_base,
                self.output_type,
                instrument,
                detector,
                list_par1,
                list_par2,
                instrument_info,
                master_table,
                **pars_cube)

#________________________________________________________________________________

            thiscube.setup_cube() # basic checks and get roi size

# find the min & max final coordinates of cube: map each slice to cube
# add any dither offsets, then find the min & max value in each dimension

            thiscube.setup_ifucube_wcs()
#________________________________________________________________________________
# build the IFU Cube

# If single = True: map each file to output grid and return single mapped file
#to output grid
# this option is used for background matching and outlier rejection
            if self.single:
                self.output_file = None
                cube_container = thiscube.build_ifucube_single()
                self.log.info("Number of Single IFUCube models returned %i ",
                              len(Final_IFUCube))

# Else standard IFU cube building
            else:
                result = thiscube.build_ifucube()
                cube_container.append(result)

            if self.debug_pixel == 1:
                self.spaxel_debug.close()
        for cube in cube_container:
            footprint = cube.meta.wcs.footprint(axis_type="spatial")
            update_s_region_keyword(cube, footprint)

        return cube_container

#********************************************************************************
class InputFileError(Exception):
    pass

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
    ValidChannel = ['1', '2', '3', '4', 'ALL']
    ValidSubChannel = ['SHORT', 'MEDIUM', 'LONG', 'ALL']
    ValidFWA = ['F070LP', 'F100LP', 'F100LP', 'F170LP',
                    'F170LP', 'F290LP', 'F290LP', 'CLEAR', 'ALL']
    ValidGWA = ['G140M', 'G140H', 'G140M', 'G140H', 'G235M', 'G235H',
                    'G395M', 'G395H', 'PRISM', 'ALL']

#________________________________________________________________________________
    # for MIRI we can set the channel
# if set to ALL then let the DetermineCubeCoverage figure out the data we have and set
# self.channel to empty
    if self.channel == 'ALL':
        self.channel = ''

    if self.channel:  # self.channel is false if it is empty
        if not self.single:
            self.output_type = 'user'
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
            else:
                raise ErrorInvalidParameter("Invalid Channel %s", ch)
# remove duplicates if needed
        self.pars_input['channel'] = list(set(self.pars_input['channel']))

#________________________________________________________________________________
    # for MIRI we can set the subchannel
# if set to ALL then let the DetermineCubeCoverage figure out the data we have and set
# self.subchannel = empty

    if self.subchannel == 'ALL':
        self.subchannel = ''

    if self.subchannel: #  not empty it has been set
        if not self.single:
            self.output_type = 'user'
        subchannellist = self.subchannel.split(',')
        user_blen = len(subchannellist)
        for j in range(user_blen):
            b = subchannellist[j]
            if user_blen > 1:
                b = b.strip('[')
                b = b.strip(']')
                b = b.strip(' ')
                b = b[1:-1]
            b = str(b)
            if b in ValidSubChannel:
                self.pars_input['subchannel'].append(b)
            else:
                raise ErrorInvalidParameter("Invalid Subchannel %s", b)
# remove duplicates if needed
        self.pars_input['subchannel'] = list(set(self.pars_input['subchannel']))
#________________________________________________________________________________
    # for NIRSPEC we can set the filter
# if set to ALL then let the DetermineCubeCoverage figure out the data we have and set
# self.filter = empty
    if self.filter == 'ALL':
        self.filter = ''
    if self.filter:
        if not self.single:
            self.output_type = 'user'
        filterlist = self.filter.split(',')
        user_flen = len(filterlist)
        for j in range(user_flen):
            f = filterlist[j]
            if user_flen > 1:
                f = f.strip('[')
                f = f.strip(']')
                f = f.strip(' ')
                f = f[1:-1]
            f = str(f)
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
        if not self.single:
            self.output_type = 'user'
        gratinglist = self.grating.split(',')
        user_glen = len(gratinglist)
        for j in range(user_glen):

            g = gratinglist[j]
            if user_glen > 1:
                g = g.strip('[')
                g = g.strip(']')
                g = g.strip(' ')
                g = g[1:-1]
            g = str(g)
            if g in ValidGWA:
                self.pars_input['grating'].append(g)
            else:
                raise ErrorInvalidParameter("Invalid Grating %s", g)
# remove duplicates if needed
        self.pars_input['grating'] = list(set(self.pars_input['grating']))
#________________________________________________________________________________
class ErrorInvalidParameter(Exception):
    pass
