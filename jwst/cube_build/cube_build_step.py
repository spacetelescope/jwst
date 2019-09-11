""" This is the main ifu spectral cube building routine.
"""
from ..stpipe import Step
from .. import datamodels
from . import cube_build
from . import ifu_cube
from . import data_types
from ..assign_wcs.util import update_s_region_keyword
from . import cube_build_wcs_util

__all__ = ["CubeBuildStep"]


class CubeBuildStep (Step):
    """CubeBuildStep: Creates a 3-D spectral cube

    Notes
    -----
    This is the controlling routine for building IFU Spectral Cubes.
    It loads and sets the various input data and parameters need by
    the cube_build_step.

    This routine does the following operations:

       1. Extracts the input parameters from the cubepars reference file and
       merges them with any user-provided values.
       2. Creates the output WCS from the input images and defines the mapping
       between all the input arrays and the output array
       3. Passes the input data to the function to map all thei input data
       to the output array.
       4. Updates the output data model with correct meta data

    """

    spec = """
         channel = option('1','2','3','4','all',default='all') # Channel
         band = option('short','medium','long','all',default='all') # Band
         grating   = option('prism','g140m','g140h','g235m','g235h',g395m','g395h','all',default='all') # Grating
         filter   = option('clear','f100lp','f070lp','f170lp','f290lp','all',default='all') # Filter
         output_type = option('band','channel','grating','multi',default='band') # Type IFUcube to create.
         scale1 = float(default=0.0) # cube sample size to use for axis 1, arc seconds
         scale2 = float(default=0.0) # cube sample size to use for axis 2, arc seconds
         scalew = float(default=0.0) # cube sample size to use for axis 3, microns
         weighting = option('msm','miripsf','area',default = 'msm') # Type of weighting function
         coord_system = option('world','alpha-beta',default='world') # Output Coordinate system. Options: world or alpha-beta
         rois = float(default=0.0) # region of interest spatial size, arc seconds
         roiw = float(default=0.0) # region of interest wavelength size, microns
         weight_power = float(default=2.0) # Weighting option to use for Modified Shepard Method
         wavemin = float(default=None)  # Minimum wavelength to be used in the IFUCube
         wavemax = float(default=None)  # Maximum wavelength to be used in the IFUCube
         single = boolean(default=false) # Internal pipeline option used by mrs_imatch & outlier detection
         xdebug = integer(default=None) # debug option, x spaxel value to report information on
         ydebug = integer(default=None) # debug option, y spaxel value to report information on
         zdebug = integer(default=None) # debug option, z spaxel value to report  information on
         skip_dqflagging = boolean(default=false) # skip setting the DQ plane of the IFU
         search_output_file = boolean(default=false)
         output_use_model = boolean(default=true) # Use filenames in the output models
       """

    reference_file_types = ['cubepar', 'resol']

# ________________________________________________________________________________
    def process(self, input):
        """This is the main routine for IFU spectral cube building.

        Parameters
        ----------
        input : list of objects or str
           list of datamodels or string name of input fits file or association.
        """

        self.log.info('Starting IFU Cube Building Step')
# ________________________________________________________________________________
# For all parameters convert to a standard format
# Report read in values to screen
# ________________________________________________________________________________
        self.subchannel = self.band
        self.suffix = 's3d'  # override suffix = cube_build

        if(not self.subchannel.islower()):
            self.subchannel = self.subchannel.lower()
        if(not self.filter.islower()):
            self.filter = self.filter.lower()
        if(not self.grating.islower()):
            self.grating = self.grating.lower()
        if(not self.coord_system.islower()):
            self.coord_system = self.coord_system.lower()
        if(not self.output_type.islower()):
            self.output_type = self.output_type.lower()
        if(not self.weighting.islower()):
            self.weighting = self.weighting.lower()

        if(self.scale1 != 0.0):
            self.log.info('Input Scale of axis 1 %f', self.scale1)
        if(self.scale2 != 0.0):
            self.log.info('Input Scale of axis 2 %f', self.scale2)
        if(self.scalew != 0.0):
            self.log.info('Input wavelength scale %f  ', self.scalew)

        if self.wavemin is not None:
            self.log.info('Setting minimum wavelength of spectral cube to: %f',
                          self.wavemin)
        if self.wavemax is not None:
            self.log.info('Setting maximum wavelength of spectral cube to: %f',
                          self.wavemax)

        if self.rois != 0.0:
            self.log.info('Input Spatial ROI size %f', self.rois)
        if self.roiw != 0.0:
            self.log.info('Input Wave ROI size %f', self.roiw)

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
        # 2. world
        self.interpolation = 'pointcloud'  # initialize

        # if the weighting is area then interpolation is area
        # valid only for single band single exposure ifucbes:
        # weighting = area or interpoltion = area

        if self.weighting == 'area':
            self.interpolation = 'area'
            self.coord_system = 'alpha-beta'

        if self.coord_system == 'alpha-beta':
            self.weighting = 'area'
            self.interpolation = 'area'

        # if interpolation is point cloud then weighting can be
        # 1. MSM: modified shepard method
        # 2. miripsf - weighting for MIRI based on PSF and LSF
        if self.coord_system == 'world':
            self.interpolation = 'pointcloud'  # can not be area

        self.log.info('Input interpolation: %s', self.interpolation)
        self.log.info('Coordinate system to use: %s', self.coord_system)
        if self.interpolation == 'pointcloud':
            self.log.info('Weighting method for point cloud: %s',
                          self.weighting)
            self.log.info('Power weighting distance : %f', self.weight_power)

        if self.single:
            self.output_type = 'single'
            self.log.info('Cube Type: Single cubes ')
            self.coord_system = 'world'
            self.interpolation = 'pointcloud'
            self.weighting = 'msm'

# ________________________________________________________________________________
# read input parameters - Channel, Band (Subchannel), Grating, Filter
# ________________________________________________________________________________
        self.pars_input = {}
# the following parameters are set either by the an input parameter
# or
# if not set on the command line then from reading in the data.
        self.pars_input['channel'] = []
        self.pars_input['subchannel'] = []

        self.pars_input['filter'] = []
        self.pars_input['grating'] = []
# read_user_input:
# see if options channel, band,grating filter are set on the command lines
# if they are then self.output_type = 'user' and fill in  par_input with values
        self.read_user_input()
# ________________________________________________________________________________
# DataTypes: Read in the input data - 4 formats are allowed:
# 1. filename
# 2. single model
# 3. ASN table
# 4. model container
# figure out what type of data we have. Fill in the input_table.input_models.
# input_table.input_models is used in the rest of IFU Cube Building
# We need to do this in cube_build_step because we need to pass the data_model
# to CRDS to figure out what type of reference files to grab (MIRI or NIRSPEC)
# if the user has provided the filename - strip out .fits and pull out the
# base name. The cube_build software will attached the needed information:
#  channel, sub-channel  grating or filter to filename
# ________________________________________________________________________________
        input_table = data_types.DataTypes(input, self.single,
                                           self.output_file,
                                           self.output_dir)

        self.input_models = input_table.input_models
        self.input_filenames = input_table.filenames
        self.output_name_base = input_table.output_name

        self.pipeline = 3
        if self.output_type == 'multi' and len(self.input_filenames) == 1:
            self.pipeline = 2

# ________________________________________________________________________________
# Read in Cube Parameter Reference file
# identify what reference file has been associated with these input
        par_filename = self.get_reference_file(self.input_models[0], 'cubepar')
# Check for a valid reference file
        if par_filename == 'N/A':
            self.log.warning('No default cube parameters reference file found')
            return
# ________________________________________________________________________________
# If miripsf weight is set then set up reference file
        resol_filename = None
        if(self.weighting == 'miripsf'):
            resol_filename = self.get_reference_file(self.input_models[0], 'resol')

            if resol_filename == 'N/A':
                self.log.warning('No spectral resolution reference file found')
                self.log.warning('Run again and turn off miripsf')
                return
# ________________________________________________________________________________
# shove the input parameters in to pars to pull out in general cube_build.py

        pars = {
            'channel': self.pars_input['channel'],
            'subchannel': self.pars_input['subchannel'],
            'grating': self.pars_input['grating'],
            'filter': self.pars_input['filter'],
            'weighting': self.weighting,
            'single': self.single,
            'output_type': self.output_type}

# shove the input parameters in to pars_cube to pull out ifu_cube.py
# these parameters are related to the building a single ifucube_model

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
            'skip_dqflagging': self.skip_dqflagging,
            'xdebug': self.xdebug,
            'ydebug': self.ydebug,
            'zdebug': self.zdebug,
            'debug_pixel': self.debug_pixel,
            'spaxel_debug': self.spaxel_debug}
# ________________________________________________________________________________
# create an instance of class CubeData

        cubeinfo = cube_build.CubeData(
            self.input_models,
            self.input_filenames,
            par_filename,
            resol_filename,
            **pars)
# ________________________________________________________________________________
# cubeinfo.setup:
# read in all the input files, information from cube_pars, read in input data
# and fill in master_table holding what files are associationed with each
# ch/sub-ch or grating/filter.
# Fill in all_channel, all_subchannel,all_filter, all_grating and instrument

        result = cubeinfo.setup()
        instrument = result['instrument']
        instrument_info = result['instrument_info']
        master_table = result['master_table']
# ________________________________________________________________________________
# How many and what type of cubes will be made.
# send self.output_type, all_channel, all_subchannel, all_grating, all_filter
# return number of cubes and for each cube the fill in
# list_pars1 (valid channel or grating) and
# list_pars2 (value subchannel or filter)

        num_cubes, cube_pars = cubeinfo.number_cubes()
        if not self.single:
            self.log.info('Number of ifucubes produced by a this run %i',
                          num_cubes)

        # ModelContainer of ifucubes
        cube_container = datamodels.ModelContainer()

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
                list_par1,
                list_par2,
                instrument_info,
                master_table,
                **pars_cube)

# ________________________________________________________________________________
            thiscube.check_ifucube()  # basic checks

# Based on channel/subchannel or grating/prism find:
# spatial spaxel size, min wave, max wave
# Set linear_wavelength to true or false depending on which type of IFUcube is
# being created.
# If linear wavelength (single band) single values are assigned to:
# rois, roiw,  weight_power, softrad
# If not linear wavelength (multi bands) an arrays are created  for:
# rois, roiw, weight_power, softrad

            thiscube.determine_cube_parameters()

            cube_footprint = thiscube.setup_ifucube_wcs()
# _______________________________________________________________________________
# build the IFU Cube

# If single = True: map each file to output grid and return single mapped file
# to output grid
# This option is used for background matching and outlier rejection
            if self.single:
                self.output_file = None
                cube_container = thiscube.build_ifucube_single()
                self.log.info("Number of Single IFUCube models returned %i ",
                              len(cube_container))

# Else standard IFU cube building
            else:
                result = thiscube.build_ifucube()
                # footprint_old = result.meta.wcs.footprint(axis_type="spatial")
                # print('old footprint',footprint_old)
                footprint_archive = cube_build_wcs_util.form_archive_footprint(cube_footprint)
                # print('new footprint', footprint_archive)
                update_s_region_keyword(result, footprint_archive)
                cube_container.append(result)
            if self.debug_pixel == 1:
                self.spaxel_debug.close()

        return cube_container
# ******************************************************************************

    def read_user_input(self):
        """Read user input options for channel, subchannel, filter, or grating

        Determine if any of the input paramters channel, band, filter or
        grating have been set. If they have been check fill in input_pars
        dictionary.

        Parameters
        ----------
        none

        Notes
        ------
        This routine updates the dictionary self.pars_input with any user
        provided inputs. In particular it sets pars_input['channel'],
        pars_input['sub_channel'], pars_input['grating'], and
        pars_input['filter'] with user provided values.
        """

        valid_channel = ['1', '2', '3', '4', 'all']
        valid_subchannel = ['short', 'medium', 'long', 'all']

        valid_fwa = ['f070lp', 'f100lp',
                     'g170lp', 'f290lp', 'clear', 'all']
        valid_gwa = ['g140m', 'g140h', 'g235m', 'g235h',
                     'g395m', 'g395h', 'prism', 'all']
# ________________________________________________________________________________
# For MIRI we can set the channel.
# If channel is  set to 'all' then let the determine_band_coverage figure out
# which channels are covered by the data.
        if self.channel == 'all':
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

                if ch in valid_channel:
                    self.pars_input['channel'].append(ch)
# remove duplicates if needed
            self.pars_input['channel'] = list(set(self.pars_input['channel']))
# ________________________________________________________________________________
# For MIRI we can set the subchannel
# if set to all then let the determine_band_coverage figure out what subchannels
# are covered by the data

        if self.subchannel == 'all':
            self.subchannel = ''

        if self.subchannel:  # not empty it has been set
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
                if b in valid_subchannel:
                    self.pars_input['subchannel'].append(b)
# remove duplicates if needed
            self.pars_input['subchannel'] = list(set(self.pars_input['subchannel']))
# ________________________________________________________________________________
# For NIRSPEC we can set the filter
# If set to all then let the determine_band_coverage figure out what filters are
# covered by the data.

        if self.filter == 'all':
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
                if f in valid_fwa:
                    self.pars_input['filter'].append(f)
# remove duplicates if needed
            self.pars_input['filter'] = list(set(self.pars_input['filter']))
# ________________________________________________________________________________
# For NIRSPEC we can set the grating
# If set to all then let the determine_band_coverage figure out what gratings are
# covered by the data
        if self.grating == 'all':
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
                if g in valid_gwa:
                    self.pars_input['grating'].append(g)
# remove duplicates if needed
            self.pars_input['grating'] = list(set(self.pars_input['grating']))
# ________________________________________________________________________________
