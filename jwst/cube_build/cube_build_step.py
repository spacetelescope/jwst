""" This is the main ifu spectral cube building routine.
"""

from jwst.datamodels import ModelContainer

from ..stpipe import Step
from . import cube_build
from . import ifu_cube
from . import data_types
from ..assign_wcs.util import update_s_region_keyword
import time

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

    class_alias = "cube_build"

    spec = """
         channel = option('1','2','3','4','all',default='all') # Channel
         band = option('short','medium','long','short-medium','short-long','medium-short', \
                'medium-long', 'long-short', 'long-medium','all',default='all') # Band
         grating   = option('prism','g140m','g140h','g235m','g235h','g395m','g395h','all',default='all') # Grating
         filter   = option('clear','f100lp','f070lp','f170lp','f290lp','all',default='all') # Filter
         output_type = option('band','channel','grating','multi',default=None) # Type IFUcube to create.
         scalexy = float(default=0.0) # cube sample size to use for axis 1 and axis2, arc seconds
         scalew = float(default=0.0) # cube sample size to use for axis 3, microns
         weighting = option('emsm','msm','drizzle',default = 'drizzle') # Type of weighting function
         coord_system = option('skyalign','world','internal_cal','ifualign',default='skyalign') # Output Coordinate system.
         rois = float(default=0.0) # region of interest spatial size, arc seconds
         roiw = float(default=0.0) # region of interest wavelength size, microns
         weight_power = float(default=2.0) # Weighting option to use for Modified Shepard Method
         wavemin = float(default=None)  # Minimum wavelength to be used in the IFUCube
         wavemax = float(default=None)  # Maximum wavelength to be used in the IFUCube
         single = boolean(default=false) # Internal pipeline option used by mrs_imatch & outlier detection
         skip_dqflagging = boolean(default=false) # skip setting the DQ plane of the IFU
         search_output_file = boolean(default=false)
         output_use_model = boolean(default=true) # Use filenames in the output models
         suffix = string(default='s3d')
         debug_spaxel = string(default='-1 -1 -1') # Default not used
       """

    reference_file_types = ['cubepar']

# ________________________________________________________________________________
    def process(self, input):
        """This is the main routine for IFU spectral cube building.

        Parameters
        ----------
        input : list of objects or str
           list of datamodels or string name of input fits file or association.
        """

        self.log.info('Starting IFU Cube Building Step')

        t0 = time.time()
# ________________________________________________________________________________
# For all parameters convert to a standard format
# Report read in values to screen
# ________________________________________________________________________________
        self.subchannel = self.band

        if not self.subchannel.islower():
            self.subchannel = self.subchannel.lower()
        if not self.filter.islower():
            self.filter = self.filter.lower()
        if not self.grating.islower():
            self.grating = self.grating.lower()
        if not self.coord_system.islower():
            self.coord_system = self.coord_system.lower()

        if not self.weighting.islower():
            self.weighting = self.weighting.lower()

        if self.scalexy != 0.0:
            self.log.info(f'Input Scale of axis 1 and 2 {self.scalexy}')
        if self.scalew != 0.0:
            self.log.info(f'Input wavelength scale {self.scalew}')

        if self.wavemin is not None:
            self.log.info(f'Setting minimum wavelength of spectral cube to: {self.wavemin}')
        if self.wavemax is not None:
            self.log.info(f'Setting maximum wavelength of spectral cube to: {self.wavemax}')

        if self.rois != 0.0:
            self.log.info(f'Input Spatial ROI size {self.rois}')
        if self.roiw != 0.0:
            self.log.info(f'Input Wave ROI size {self.roiw}')

        # valid coord_system:
        # 1. skyalign (ra dec) (aka world)
        # 2. ifualign (ifu cube aligned with slicer plane/ MRS local coord system)
        # 3. internal_cal (local IFU - ifu cubes built in local IFU system)
        if self.coord_system == 'world':
            self.coord_system = 'skyalign'

        self.interpolation = 'pointcloud'  # initialize

        # coord system = internal_cal only option for weighting = area ONLY USED FOR NIRSPEC
        if self.coord_system == 'internal_cal':
            self.interpolation = 'area'
            self.weighting = 'emsm'

        # if interpolation is point cloud then weighting can be
        # 1. MSM: modified Shepard method
        # 2. EMSM

        if self.coord_system == 'skyalign':
            self.interpolation = 'pointcloud'

        if self.coord_system == 'ifualign':
            self.interpolation = 'pointcloud'

        if self.weighting == 'drizzle':
            self.interpolation = 'drizzle'

        self.log.info(f'Input interpolation: {self.interpolation}')
        self.log.info(f'Coordinate system to use: {self.coord_system}')
        if self.interpolation == 'pointcloud':
            self.log.info(f'Weighting method for point cloud: {self.weighting}')
            if self.weight_power != 0:
                self.log.info(f'Power weighting distance: {self.weight_power}')

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

        # including values in pars_input that could get updated in cube_build_step.py

        self.pars_input['coord_system'] = self.coord_system

        if self.single:
            self.pars_input['output_type'] = 'single'
            self.log.info('Cube Type: Single cubes')
            self.pars_input['coord_system'] = 'skyalign'

            # Don't allow anything but drizzle, msm, or emsm weightings
            if self.weighting not in ['msm', 'emsm', 'drizzle']:
                self.weighting = 'drizzle'

            if self.weighting == 'drizzle':
                self.interpolation = 'drizzle'

            if self.weighting == 'msm':
                self.interpolation = 'pointcloud'

            if self.weighting == 'emsm':
                self.interpolation = 'pointcloud'

        # read_user_input:
        # see if options channel, band,grating filter are set on the command lines
        # if they are then self.pars_input['output_type'] = 'user' and fill in  par_input with values
        self.read_user_input()
# ________________________________________________________________________________
# DataTypes
# Read in the input data - 4 formats are allowed:
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
        self.output_name_base = input_table.output_name
# ________________________________________________________________________________

# Read in the first input model to determine with instrument we have
# output type is by default 'Channel' for MIRI and 'Band' for NIRSpec
        instrument = self.input_models[0].meta.instrument.name.upper()
        if self.output_type is None:
            if instrument == 'NIRSPEC':
                self.output_type = 'band'

            elif instrument == 'MIRI':
                self.output_type = 'channel'
        self.pars_input['output_type'] = self.output_type
        
        self.log.info(f'Setting output type to: {self.output_type}')
         
# Read in Cube Parameter Reference file
# identify what reference file has been associated with these input

        par_filename = self.get_reference_file(self.input_models[0], 'cubepar')
        # Check for a valid reference file
        if par_filename == 'N/A':
            self.log.warning('No default cube parameters reference file found')
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
            'output_type': self.pars_input['output_type']}

        # shove the input parameters in to pars_cube to pull out ifu_cube.py
        # these parameters are related to the building a single ifucube_model

        pars_cube = {
            'scalexy': self.scalexy,
            'scalew': self.scalew,
            'interpolation': self.interpolation,
            'weighting': self.weighting,
            'weight_power': self.weight_power,
            'coord_system': self.pars_input['coord_system'],
            'rois': self.rois,
            'roiw': self.roiw,
            'wavemin': self.wavemin,
            'wavemax': self.wavemax,
            'skip_dqflagging': self.skip_dqflagging,
            'suffix': self.suffix,
            'debug_spaxel': self.debug_spaxel}

# ________________________________________________________________________________
# create an instance of class CubeData

        cubeinfo = cube_build.CubeData(
            self.input_models,
            par_filename,
            **pars)
# ________________________________________________________________________________
# cubeinfo.setup:
# read in all the input files, information from cube_pars, read in input data
# and fill in master_table holding what files are associated with each
# ch/sub-ch or grating/filter.
# Fill in all_channel, all_subchannel,all_filter, all_grating and instrument

        result = cubeinfo.setup()
        instrument = result['instrument']
        instrument_info = result['instrument_info']
        master_table = result['master_table']

        if instrument == 'MIRI' and self.coord_system == 'internal_cal':
            self.log.warning('The output coordinate system of internal_cal is not valid for MIRI')
            self.log.warning('use output_coord = ifualign instead')
            return
        filenames = master_table.FileMap['filename']

        self.pipeline = 3
        if self.output_type == 'multi' and len(filenames) == 1:
            self.pipeline = 2
# ________________________________________________________________________________
# How many and what type of cubes will be made.
# send self.pars_input['output_type'], all_channel, all_subchannel, all_grating, all_filter
# return number of cubes and for each cube the fill in
# list_pars1 (valid channel or grating) and
# list_pars2 (value subchannel or filter)

        num_cubes, cube_pars = cubeinfo.number_cubes()
        if not self.single:
            self.log.info(f'Number of IFU cubes produced by this run = {num_cubes}')

        # ModelContainer of ifucubes
        cube_container = ModelContainer()
        status_cube = 0

        # for single type cubes num_cubes always = 1, Looping over
        # bands is done in outlier detection.
        for i in range(num_cubes):
            icube = str(i + 1)
            list_par1 = cube_pars[icube]['par1']
            list_par2 = cube_pars[icube]['par2']
            thiscube = ifu_cube.IFUCubeData(
                self.pipeline,
                self.input_models,
                self.output_name_base,
                self.pars_input['output_type'],
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

            if self.coord_system == 'internal_cal':
                thiscube.determine_cube_parameters_internal()
            else:
                thiscube.determine_cube_parameters()
            thiscube.setup_ifucube_wcs()
# _______________________________________________________________________________
# build the IFU Cube

            status = 0

# If single = True: map each file to output grid and return single mapped file
# to output grid. # This option is used for background matching and outlier rejection

            if self.single:
                self.output_file = None
                cube_container = thiscube.build_ifucube_single()
                self.log.info("Number of Single IFUCube models returned %i ",
                              len(cube_container))

# Else standard IFU cube building - the result returned from build_ifucube will be 1 IFU CUBR
            else:
                result, status = thiscube.build_ifucube()

                # check if cube_build failed
                # **************************
                if status == 1:
                    status_cube = 1

                cube_container.append(result)
                del result
            del thiscube

        # irrelevant WCS keywords we will remove from final product
        rm_keys = ['v2_ref', 'v3_ref', 'ra_ref', 'dec_ref', 'roll_ref',
                   'v3yangle', 'vparity']

        for cube in cube_container:
            footprint = cube.meta.wcs.footprint(axis_type="spatial")
            update_s_region_keyword(cube, footprint)

            # remove certain WCS keywords that are irrelevant after combine data into IFUCubes
            for key in rm_keys:
                if key in cube.meta.wcsinfo.instance:
                    del cube.meta.wcsinfo.instance[key]
        if status_cube == 1:
            self.skip = True

        t1 = time.time()
        self.log.debug(f'Time to build all cubes {t1-t0}')

        if status_cube == 1:
            self.skip = True

        return cube_container
# ******************************************************************************

    def read_user_input(self):
        """Read user input options for channel, subchannel, filter, or grating"""

        # Determine if any of the input parameters channel, band, filter or
        # grating have been set.

        # This routine updates the dictionary self.pars_input with any user
        # provided inputs. In particular it sets pars_input['channel'],
        # pars_input['sub_channel'], pars_input['grating'], and
        # pars_input['filter'] with user provided values.

        valid_channel = ['1', '2', '3', '4', 'all']
        valid_subchannel = ['short', 'medium', 'long', 'all', 'short-medium', 'short-long',
                            'medium-short', 'medium-long', 'long-short', 'long-medium']

        valid_fwa = ['f070lp', 'f100lp',
                     'g170lp', 'f290lp', 'clear', 'opaque', 'all']
        valid_gwa = ['g140m', 'g140h', 'g235m', 'g235h',
                     'g395m', 'g395h', 'prism', 'all']
# ________________________________________________________________________________
# For MIRI we can set the channel.
# If channel is  set to 'all' then let the determine_band_coverage figure out
# which channels are covered by the data.

        if self.channel == 'all':
            self.pars_input['channel'].append('all')
        else:  # user has set value
            if not self.single:
                self.pars_input['output_type'] = 'user'
            channellist = self.channel.split(',')
            user_clen = len(channellist)
            for j in range(user_clen):
                ch = channellist[j]
                if user_clen > 1:
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
            self.pars_input['subchannel'].append('all')
        else:  # user has set value
            if not self.single:
                self.pars_input['output_type'] = 'user'
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
            self.pars_input['filter'].append('all')
        else:   # User has set value
            if not self.single:
                self.pars_input['output_type'] = 'user'
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
            self.pars_input['grating'].append('all')
        else:    # user has set value
            if not self.single:
                self.pars_input['output_type'] = 'user'
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
