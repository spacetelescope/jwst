# Routines used for building cubes
from .. import datamodels
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

##################################################################################
class FileTable():
    """
    Dictionary that maps the input files to correct:
    MIRI - channel/subchannel
    NIRSPEC - grating/filter

    Parameters
    ----------

    """
    def __init__(self):

        self.FileMap = {}
        self.FileMap['MIRI'] = {}

        self.FileMap['MIRI']['1'] = {}
        self.FileMap['MIRI']['1']['short'] = []
        self.FileMap['MIRI']['1']['medium'] = []
        self.FileMap['MIRI']['1']['long'] = []

        self.FileMap['MIRI']['2'] = {}
        self.FileMap['MIRI']['2']['short'] = []
        self.FileMap['MIRI']['2']['medium'] = []
        self.FileMap['MIRI']['2']['long'] = []

        self.FileMap['MIRI']['3'] = {}
        self.FileMap['MIRI']['3']['short'] = []
        self.FileMap['MIRI']['3']['medium'] = []
        self.FileMap['MIRI']['3']['long'] = []

        self.FileMap['MIRI']['4'] = {}
        self.FileMap['MIRI']['4']['short'] = []
        self.FileMap['MIRI']['4']['medium'] = []
        self.FileMap['MIRI']['4']['long'] = []

        self.FileMap['NIRSPEC'] = {}
        self.FileMap['NIRSPEC']['prism'] = {}
        self.FileMap['NIRSPEC']['prism']['clear'] = []

        self.FileMap['NIRSPEC']['g140m'] = {}
        self.FileMap['NIRSPEC']['g140m']['f070lp'] = []
        self.FileMap['NIRSPEC']['g140m']['f100lp'] = []

        self.FileMap['NIRSPEC']['g140h'] = {}
        self.FileMap['NIRSPEC']['g140h']['f070lp'] = []
        self.FileMap['NIRSPEC']['g140h']['f100lp'] = []

        self.FileMap['NIRSPEC']['g235m'] = {}
        self.FileMap['NIRSPEC']['g235m']['f170lp'] = []

        self.FileMap['NIRSPEC']['g235h'] = {}
        self.FileMap['NIRSPEC']['g235h']['f170lp'] = []

        self.FileMap['NIRSPEC']['g395m'] = {}
        self.FileMap['NIRSPEC']['g395m']['f290lp'] = []

        self.FileMap['NIRSPEC']['g395h'] = {}
        self.FileMap['NIRSPEC']['g395h']['f290lp'] = []

        self.FileOffset = {}
        self.FileOffset['1'] = {}
        self.FileOffset['1']['short'] = {}
        self.FileOffset['1']['short']['C1'] = []
        self.FileOffset['1']['short']['C2'] = []
        self.FileOffset['1']['medium'] = {}
        self.FileOffset['1']['medium']['C1'] = []
        self.FileOffset['1']['medium']['C2'] = []
        self.FileOffset['1']['long'] = {}
        self.FileOffset['1']['long']['C1'] = []
        self.FileOffset['1']['long']['C2'] = []

        self.FileOffset['2'] = {}
        self.FileOffset['2']['short'] = {}
        self.FileOffset['2']['short']['C1'] = []
        self.FileOffset['2']['short']['C2'] = []
        self.FileOffset['2']['medium'] = {}
        self.FileOffset['2']['medium']['C1'] = []
        self.FileOffset['2']['medium']['C2'] = []
        self.FileOffset['2']['long'] = {}
        self.FileOffset['2']['long']['C1'] = []
        self.FileOffset['2']['long']['C2'] = []

        self.FileOffset['3'] = {}
        self.FileOffset['3']['short'] = {}
        self.FileOffset['3']['short']['C1'] = []
        self.FileOffset['3']['short']['C2'] = []
        self.FileOffset['3']['medium'] = {}
        self.FileOffset['3']['medium']['C1'] = []
        self.FileOffset['3']['medium']['C2'] = []
        self.FileOffset['3']['long'] = {}
        self.FileOffset['3']['long']['C1'] = []
        self.FileOffset['3']['long']['C2'] = []

        self.FileOffset['4'] = {}
        self.FileOffset['4']['short'] = {}
        self.FileOffset['4']['short']['C1'] = []
        self.FileOffset['4']['short']['C2'] = []
        self.FileOffset['4']['medium'] = {}
        self.FileOffset['4']['medium']['C1'] = []
        self.FileOffset['4']['medium']['C2'] = []
        self.FileOffset['4']['long'] = {}
        self.FileOffset['4']['long']['C1'] = []
        self.FileOffset['4']['long']['C2'] = []

        self.FileOffset['prism'] = {}
        self.FileOffset['prism']['clear'] = {}
        self.FileOffset['prism']['clear']['C1'] = []
        self.FileOffset['prism']['clear']['C2'] = []

        self.FileOffset['g140m'] = {}
        self.FileOffset['g140m']['f070lp'] = {}
        self.FileOffset['g140m']['f070lp']['C1'] = []
        self.FileOffset['g140m']['f070lp']['C2'] = []
        self.FileOffset['g140m']['f100lp'] = {}
        self.FileOffset['g140m']['f100lp']['C1'] = []
        self.FileOffset['g140m']['f100lp']['C2'] = []

        self.FileOffset['g140h'] = {}
        self.FileOffset['g140h']['f070lp'] = {}
        self.FileOffset['g140h']['f070lp']['C1'] = []
        self.FileOffset['g140h']['f070lp']['C2'] = []
        self.FileOffset['g140h']['f100lp'] = {}
        self.FileOffset['g140h']['f100lp']['C1'] = []
        self.FileOffset['g140h']['f100lp']['C2'] = []

        self.FileOffset['g235m'] = {}
        self.FileOffset['g235m']['f170lp'] = {}
        self.FileOffset['g235m']['f170lp']['C1'] = []
        self.FileOffset['g235m']['f170lp']['C2'] = []

        self.FileOffset['g235h'] = {}
        self.FileOffset['g235h']['f170lp'] = {}
        self.FileOffset['g235h']['f170lp']['C1'] = []
        self.FileOffset['g235h']['f170lp']['C2'] = []

        self.FileOffset['g395m'] = {}
        self.FileOffset['g395m']['f290lp'] = {}
        self.FileOffset['g395m']['f290lp']['C1'] = []
        self.FileOffset['g395m']['f290lp']['C2'] = []

        self.FileOffset['g395h'] = {}
        self.FileOffset['g395h']['f290lp'] = {}
        self.FileOffset['g395h']['f290lp']['C1'] = []
        self.FileOffset['g395h']['f290lp']['C2'] = []


#********************************************************************************
    def set_file_table(self,
                       input_models,
                       input_filenames,
                       input_ra_offset,
                       input_dec_offset):
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


        Returns
        -------
        MasterTable filled in with files needed
        num: number of files to create cube from
        detector

        """
        num = 0
        num = len(input_filenames)
#________________________________________________________________________________
# Loop over input list of files and assign fill in the MasterTable with filename
# for the correct (channel-subchannel) or (grating-subchannel)
        for i in range(num):

            ifile = input_filenames[i]
            input = input_models[i]

        # Open the input data model & Fill in the FileMap information

            with datamodels.IFUImageModel(input) as input_model:

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
                    subchannel = input_model.meta.instrument.band.lower()

            #________________________________________________________________________________
                    clenf = len(channel)
                    for k in range(clenf):
                        self.FileMap['MIRI'][channel[k]][subchannel].append(input_model)
                        ioffset = len(input_ra_offset)
                        if (ioffset > 0):
                            ra_offset = input_ra_offset[i]
                            dec_offset = input_dec_offset[i]
                            self.FileOffset[channel[k]][subchannel]['C1'].append(ra_offset)
                            self.FileOffset[channel[k]][subchannel]['C2'].append(dec_offset)
            #________________________________________________________________________________
                elif instrument == 'NIRSPEC':
                    fwa = input_model.meta.instrument.filter.lower()
                    gwa = input_model.meta.instrument.grating.lower()

                    self.FileMap['NIRSPEC'][gwa][fwa].append(input_model)
                else:

                    log.info('Instrument not valid for cube')

        return instrument, detector

class ErrorNoAssignWCS(Exception):
    pass
