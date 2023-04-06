"""  Basic routines used to set up IFU cubes
"""
import logging
from . import cube_build_io_util
from . import file_table
from . import instrument_defaults

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class CubeData():
    """ class CubeData holds top level information on the ifucube

    """

    def __init__(self,
                 input_models,
                 par_filename,
                 **pars):
        """ Initialize the high level of information for the ifu cube

        The class CubeData holds information on what type of cube is being
        created and how the cube is to be constructed. This information includes:
        1. Bands covered by data (channel/sub-channel or filter/grating)
        2. If the IFU cube is a single or multiple band
        3. What instrument the data is for
        4. weighting function to use to construct the ifu cube

        Parameters
        ----------
        input_models : list of data models
        par_filename: str
          cube parameter reference filename
        pars : dictionary holding top level cube parameters
        """

        self.input_models = input_models
        self.par_filename = par_filename
        self.single = pars.get('single')
        self.channel = pars.get('channel')
        self.subchannel = pars.get('subchannel')
        self.grating = pars.get('grating')
        self.filter = pars.get('filter')
        self.weighting = pars.get('weighting')
        self.output_type = pars.get('output_type')
        self.instrument = None

        self.all_channel = []
        self.all_subchannel = []
        self.all_grating = []
        self.all_filter = []

        self.output_name = ''
# _____________________________________________________________________________

    def setup(self):
        """ Set up IFU Cube. Determine band coverage and read in reference files

        Read in the input_models and fill in the dictionary, master_table, that
        stores the data  for each channel/subchannel or grating/filter.

        If the channel/subchannel or grating/filter are not set by the user,
        then determine which ones are found in the data

        Read in necessary reference data:
        * cube parameter reference file

        This routine fills in the instrument_info dictionary, which holds the
        default spatial and spectral size of the output cube, as well as,
        the region of influence size in the spatial and spectral dimensions.
        The dictionary master_table is also filled in. This dictionary contains
        list of datamodels for each band (channel/subchannel or grating/filter
        combination).
        """
# _____________________________________________________________________________
# Read in the input data (association table or single file)
# Fill in MasterTable   based on Channel/Subchannel  or filter/grating
# ______________________________________________________________________________
        master_table = file_table.FileTable()
        instrument = master_table.set_file_table(self.input_models)
# _______________________________________________________________________________
# find out how many files are in the association table or if it is an single
# file store the input_filenames and input_models

        self.instrument = instrument
# _______________________________________________________________________________
# Determine which channels/subchannels or filter/grating will be covered by the
# spectral cube.
# fills in band_channel, band_subchannel, band_grating, band_filer
# _______________________________________________________________________________
        self.determine_band_coverage(master_table)
# _______________________________________________________________________________
# instrument_defaults is an  dictionary class that holds default parameters for
# each band for the different instruments
# _______________________________________________________________________________
        instrument_info = instrument_defaults.InstrumentInfo()
# -------------------------------------------------------------------------------
# Read the cube pars reference file
        log.info('Reading cube parameter file %s', self.par_filename)
        cube_build_io_util.read_cubepars(self.par_filename,
                                         self.instrument,
                                         self.weighting,
                                         self.all_channel,
                                         self.all_subchannel,
                                         self.all_grating,
                                         self.all_filter,
                                         instrument_info)
# -------------------------------------------------------------------------------
        self.instrument_info = instrument_info
# _______________________________________________________________________________
# Set up values to return and access for other parts of cube_build

        self.master_table = master_table

        return {'instrument': self.instrument,
                'instrument_info': self.instrument_info,
                'master_table': self.master_table}

# ********************************************************************************

    def determine_band_coverage(self, master_table):
        """Either by user parameters or the data itself, determine which bands
        cover by cube

        To determine which bands are covered by the output cube:
        1. If MIRI data then check if the user has set either the channel or
        band to use. If these have not been set then read in the data to find
        out which bands the input data  covered.
        2. If NIRSPEC data then check if the user has set the grating and
        filter to use. If these have not been set then read in the data to find
        out which bands the input data covers.

        Function to determine which files contain channels and subchannels
        are used in the creation of the cubes.
        For MIRI The channels  to be used are set by the association and the
        subchannels are  determined from the data.

        Parameters
        ----------
        master_table: dictionary
          A dictionary for each band that contains of list of datamodels
          covering that band.

        Raises
        ------
        ErrorNoChannels
           The user selected channels are not in the data
        ErrorNoSubChannels
           The user selected subchannels are not in the data
        ErrorNoGratings
           The user selected gratings are not in the data
        ErrorNoFilters
           The user selected filters are not in the data
        ErrorMissingParameter
           The user selected grating but not filter or vice versa
        """
# ________________________________________________________________________________
# IF INSTRUMENT = MIRI
# loop over the file names
        if self.instrument == 'MIRI':
            valid_channel = ['1', '2', '3', '4']
            valid_subchannel = ['short', 'medium', 'long', 'short-medium', 'short-long',
                                'medium-short', 'medium-long', 'long-short', 'long-medium']

            nchannels = len(valid_channel)
            nsubchannels = len(valid_subchannel)
            # _______________________________________________________________________________
            # for MIRI we can set the channel and subchannel
            # check if default of 'all' should be used or if the user has set option
            # if default 'all' is used then we read in the data to figure which channels &
            # sub channels we need to use

            user_clen = len(self.channel)
            user_slen = len(self.subchannel)

            if user_clen == 1 and self.channel[0] == 'all':
                user_clen = 0

            if user_slen == 1 and self.subchannel[0] == 'all':
                user_slen = 0

            # _______________________________________________________________________________
            for i in range(nchannels):
                for j in range(nsubchannels):
                    nfiles = len(master_table.FileMap['MIRI'][valid_channel[i]][valid_subchannel[j]])
                    if nfiles > 0:
                        # ______________________________________________________
                        # neither parameters are set
                        if user_clen == 0 and user_slen == 0:
                            self.all_channel.append(valid_channel[i])
                            self.all_subchannel.append(valid_subchannel[j])
                        # channel was set by user but not sub-channel
                        elif user_clen != 0 and user_slen == 0:
                            if valid_channel[i] in self.channel:
                                self.all_channel.append(valid_channel[i])
                                self.all_subchannel.append(valid_subchannel[j])
                        # sub-channel was set by user but not channel
                        elif user_clen == 0 and user_slen != 0:
                            if valid_subchannel[j] in self.subchannel:
                                self.all_channel.append(valid_channel[i])
                                self.all_subchannel.append(valid_subchannel[j])
                        # both parameters set
                        else:
                            if (valid_channel[i] in self.channel and
                                    valid_subchannel[j] in self.subchannel):
                                self.all_channel.append(valid_channel[i])
                                self.all_subchannel.append(valid_subchannel[j])

            log.info('The desired cubes cover the MIRI Channels: %s',
                     self.all_channel)
            log.info('The desired cubes cover the MIRI subchannels: %s',
                     self.all_subchannel)

            number_channels = len(self.all_channel)
            number_subchannels = len(self.all_subchannel)

            if number_channels == 0:
                raise ErrorNoChannels(
                    "The cube  does not cover any channels, change channel parameter")
            if number_subchannels == 0:
                raise ErrorNoSubchannels(
                    "The cube does not cover any subchannels, change band parameter")
        # _______________________________________________________________
        if self.instrument == 'NIRSPEC':
            # 1 to 1 mapping valid_gwa[i] -> valid_fwa[i]
            valid_gwa = ['g140m', 'g140h', 'g140m', 'g140h', 'g235m',
                         'g235h', 'g395m', 'g395h', 'prism',
                         'prism', 'g140m', 'g140h', 'g235m', 'g235h', 'g395m', 'g395h']
            valid_fwa = ['f070lp', 'f070lp', 'f100lp', 'f100lp', 'f170lp',
                         'f170lp', 'f290lp', 'f290lp', 'clear',
                         'opaque', 'opaque', 'opaque', 'opaque', 'opaque', 'opaque', 'opaque']

            nbands = len(valid_fwa)

            user_glen = len(self.grating)
            user_flen = len(self.filter)
            if user_glen == 1 and self.grating[0] == 'all':
                user_glen = 0

            if user_flen == 1 and self.filter[0] == 'all':
                user_flen = 0

            # check if input filter or grating has been set
            if user_glen == 0 and user_flen != 0:
                raise ErrorMissingParameter("Filter specified, but Grating was not")

            for i in range(nbands):
                nfiles = len(master_table.FileMap['NIRSPEC'][valid_gwa[i]][valid_fwa[i]])
                if nfiles > 0:
                    # _________________________________________________
                    # neither parameters are set
                    if user_glen == 0 and user_flen == 0:
                        self.all_grating.append(valid_gwa[i])
                        self.all_filter.append(valid_fwa[i])
                    # __________________________________________________
                    # grating was set by user but not filter
                    elif user_glen != 0 and user_flen == 0:
                        if valid_gwa[i] in self.grating:
                            self.all_grating.append(valid_gwa[i])
                            self.all_filter.append(valid_fwa[i])
                    # __________________________________________________
                    # both parameters set
                    else:
                        if (valid_fwa[i] in self.filter and
                                valid_gwa[i] in self.grating):
                            self.all_grating.append(valid_gwa[i])
                            self.all_filter.append(valid_fwa[i])

            number_filters = len(self.all_filter)
            number_gratings = len(self.all_grating)

            if number_filters == 0:
                raise ErrorNoFilters("The cube does not cover any filters")
            if number_gratings == 0:
                raise ErrorNoGratings("The cube does not cover any gratings")
# ______________________________________________________________________

    def number_cubes(self):
        """Determine the number of IFUcubes to created based on:
        Type of cube (single band, multiple bands, or Single mode)
        """
        num_cubes = 0
        cube_pars = {}
# ______________________________________________________________________
# MIRI
# ______________________________________________________________________
        if self.instrument == 'MIRI':
            band_channel = self.all_channel
            band_subchannel = self.all_subchannel

# user, single, or multi
            if (self.output_type == 'user' or self.output_type == 'single' or
                    self.output_type == 'multi'):

                if self.output_type == 'multi':
                    log.info('Output IFUcube are constructed from all the data ')
                if self.single:
                    log.info('Single = true, creating a set of single exposures mapped' +
                             ' to output IFUCube coordinate system')
                if self.output_type == 'user':
                    log.info('The user has selected the type of IFU cube to make')

                num_cubes = 1
                cube_pars['1'] = {}
                cube_pars['1']['par1'] = {}
                cube_pars['1']['par2'] = {}
                cube_pars['1']['par1'] = self.all_channel
                cube_pars['1']['par2'] = self.all_subchannel

# default band cubes
            if self.output_type == 'band':
                log.info('Output Cubes are single channel, single sub-channel IFU Cubes')

                for i in range(len(band_channel)):
                    num_cubes = num_cubes + 1
                    cube_no = str(num_cubes)
                    cube_pars[cube_no] = {}
                    cube_pars[cube_no]['pars1'] = {}
                    cube_pars[cube_no]['pars2'] = {}
                    this_channel = []
                    this_subchannel = []
                    this_channel.append(band_channel[i])
                    this_subchannel.append(band_subchannel[i])
                    cube_pars[cube_no]['par1'] = this_channel
                    cube_pars[cube_no]['par2'] = this_subchannel

# default channel cubes
            if self.output_type == 'channel':
                log.info('Output cubes are single channel and all subchannels in data')
                num_cubes = 0
                channel_no_repeat = list(set(band_channel))

                for i in channel_no_repeat:
                    num_cubes = num_cubes + 1
                    cube_no = str(num_cubes)
                    cube_pars[cube_no] = {}
                    cube_pars[cube_no]['pars1'] = {}
                    cube_pars[cube_no]['pars2'] = {}
                    this_channel = []
                    this_subchannel = []
                    for k, j in enumerate(band_channel):
                        if j == i:
                            this_subchannel.append(band_subchannel[k])
                            this_channel.append(i)
                    cube_pars[cube_no]['par1'] = this_channel
                    cube_pars[cube_no]['par2'] = this_subchannel

# ______________________________________________________________________
# NIRSPEC
# ______________________________________________________________________
        if self.instrument == 'NIRSPEC':

            band_grating = self.all_grating
            band_filter = self.all_filter

            if (self.output_type == 'user' or self.output_type == 'single' or
                    self.output_type == 'multi'):
                if self.output_type == 'multi':
                    log.info('Output IFUcube are constructed from all the data ')
                if self.single:
                    log.info('Single = true, creating a set of single exposures' +
                             ' mapped to output IFUCube coordinate system')
                if self.output_type == 'user':
                    log.info('The user has selected the type of IFU cube to make')

                num_cubes = 1
                cube_pars['1'] = {}
                cube_pars['1']['par1'] = {}
                cube_pars['1']['par2'] = {}
                cube_pars['1']['par1'] = self.all_grating
                cube_pars['1']['par2'] = self.all_filter

# default band cubes
            if self.output_type == 'band':
                log.info('Output Cubes are single grating, single filter IFU Cubes')
                for i in range(len(band_grating)):
                    num_cubes = num_cubes + 1
                    cube_no = str(num_cubes)
                    cube_pars[cube_no] = {}
                    cube_pars[cube_no]['pars1'] = {}
                    cube_pars[cube_no]['pars2'] = {}
                    this_grating = []
                    this_filter = []
                    this_grating.append(band_grating[i])
                    this_filter.append(band_filter[i])
                    cube_pars[cube_no]['par1'] = this_grating
                    cube_pars[cube_no]['par2'] = this_filter
# default grating cubes
            if self.output_type == 'grating':
                log.info('Output cubes are single grating & all filters in data')
                num_cubes = 0
                for i in range(len(band_grating)):
                    num_cubes = num_cubes + 1
                    cube_no = str(num_cubes)
                    cube_pars[cube_no] = {}
                    cube_pars[cube_no]['pars1'] = {}
                    cube_pars[cube_no]['pars2'] = {}
                    this_grating = []
                    this_filter = band_subchannel
                    this_grating.append(i)
                    cube_pars[cube_no]['par1'] = this_grating
                    cube_pars[cube_no]['par2'] = this_filter

        self.num_cubes = num_cubes
        self.cube_pars = cube_pars
        return self.num_cubes, self.cube_pars
# _____________________________________________________________________________


class ErrorNoChannels(Exception):
    """ Raises Exception if the user selected channels are not in the data
    """
    pass


class ErrorNoSubchannels(Exception):
    """ Raises Exception if the user selected subchannels are not in the data
    """
    pass


class ErrorNoFilters(Exception):
    """ Raises Exception if the user selected filters are not in the data
    """
    pass


class ErrorNoGratings(Exception):
    """ Raises Exception if the user selected gratings are not in the data
    """
    pass


class ErrorMissingParameter(Exception):
    """ Raises Exception if provided grating but not filter or vice versa
    """
    pass
