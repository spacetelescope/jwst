""" Dictionary holding defaults for cube_build
"""
from .. import datamodels
from .. assign_wcs.util import  wcs_bbox_from_shape
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class FileTable():
    """ Dictionary contains defaults for MIRI and NIRSPEC data
    """
    def __init__(self):

        self.FileMap = {}
        self.FileMap['MIRI'] = {}

#        self.FileMap['MIRI']['1'] = {}
#        self.FileMap['MIRI']['1']['short'] = []
#        self.FileMap['MIRI']['1']['medium'] = []
#        self.FileMap['MIRI']['1']['long'] = []
        self.FileMap['MIRI']['1'] = {}
        self.FileMap['MIRI']['1']['short'] = {}
        self.FileMap['MIRI']['1']['short']['footprint'] = []
        self.FileMap['MIRI']['1']['short']['file'] = []
        self.FileMap['MIRI']['1']['medium'] = {}
        self.FileMap['MIRI']['1']['medium']['footprint'] = []
        self.FileMap['MIRI']['1']['medium']['file'] = []
        self.FileMap['MIRI']['1']['long'] = {}
        self.FileMap['MIRI']['1']['long']['footprint'] = []
        self.FileMap['MIRI']['1']['long']['file'] = []

        self.FileMap['MIRI']['2'] = {}
        self.FileMap['MIRI']['2']['short'] = {}
        self.FileMap['MIRI']['2']['short']['footprint'] = []
        self.FileMap['MIRI']['2']['short']['file'] = []
        self.FileMap['MIRI']['2']['medium'] = {}
        self.FileMap['MIRI']['2']['medium']['footprint'] = []
        self.FileMap['MIRI']['2']['medium']['file'] = []
        self.FileMap['MIRI']['2']['long'] = {}
        self.FileMap['MIRI']['2']['long']['footprint'] = []
        self.FileMap['MIRI']['2']['long']['file'] = []

        self.FileMap['MIRI']['3'] = {}
        self.FileMap['MIRI']['3']['short'] = {}
        self.FileMap['MIRI']['3']['short']['footprint'] = []
        self.FileMap['MIRI']['3']['short']['file'] = []
        self.FileMap['MIRI']['3']['medium'] = {}
        self.FileMap['MIRI']['3']['medium']['footprint'] = []
        self.FileMap['MIRI']['3']['medium']['file'] = []
        self.FileMap['MIRI']['3']['long'] = {}
        self.FileMap['MIRI']['3']['long']['footprint'] = []
        self.FileMap['MIRI']['3']['long']['file'] = []

        self.FileMap['MIRI']['4'] = {}
        self.FileMap['MIRI']['4']['short'] = {}
        self.FileMap['MIRI']['4']['short']['footprint'] = []
        self.FileMap['MIRI']['4']['short']['file'] = []
        self.FileMap['MIRI']['4']['medium'] = {}
        self.FileMap['MIRI']['4']['medium']['footprint'] = []
        self.FileMap['MIRI']['4']['medium']['file'] = []
        self.FileMap['MIRI']['4']['long'] = {}
        self.FileMap['MIRI']['4']['long']['footprint'] = []
        self.FileMap['MIRI']['4']['long']['file'] = []

        self.FileMap['NIRSPEC'] = {}
        self.FileMap['NIRSPEC']['prism'] = {}
        self.FileMap['NIRSPEC']['prism']['clear'] = {}
        self.FileMap['NIRSPEC']['prism']['clear']['footprint'] = []
        self.FileMap['NIRSPEC']['prism']['clear']['file'] = []

        self.FileMap['NIRSPEC']['g140m'] = {}
        self.FileMap['NIRSPEC']['g140m']['f070lp'] = {}
        self.FileMap['NIRSPEC']['g140m']['f070lp']['footprint'] = []
        self.FileMap['NIRSPEC']['g140m']['f070lp']['file'] = []

#        self.FileMap['NIRSPEC']['g140m']['f100lp'] = []
        self.FileMap['NIRSPEC']['g140m']['f100lp'] = {}
        self.FileMap['NIRSPEC']['g140m']['f100lp']['file'] = []
        self.FileMap['NIRSPEC']['g140m']['f100lp']['footprint'] = []

        self.FileMap['NIRSPEC']['g140h'] = {}
        self.FileMap['NIRSPEC']['g140h']['f070lp'] = {}
        self.FileMap['NIRSPEC']['g140h']['f070lp']['footprint'] = []
        self.FileMap['NIRSPEC']['g140h']['f070lp']['file'] = []

        self.FileMap['NIRSPEC']['g140h']['f100lp'] = {}
        self.FileMap['NIRSPEC']['g140h']['f100lp']['footprint'] = []
        self.FileMap['NIRSPEC']['g140h']['f100lp']['file'] = []

        self.FileMap['NIRSPEC']['g235m'] = {}
        self.FileMap['NIRSPEC']['g235m']['f170lp'] = {}
        self.FileMap['NIRSPEC']['g235m']['f170lp']['footprint'] = []
        self.FileMap['NIRSPEC']['g235m']['f170lp']['file'] = []

        self.FileMap['NIRSPEC']['g235h'] = {}
        self.FileMap['NIRSPEC']['g235h']['f170lp'] = {}
        self.FileMap['NIRSPEC']['g235h']['f170lp']['footprint'] = []
        self.FileMap['NIRSPEC']['g235h']['f170lp']['file'] = []

        self.FileMap['NIRSPEC']['g395m'] = {}
        self.FileMap['NIRSPEC']['g395m']['f290lp'] = {}
        self.FileMap['NIRSPEC']['g395m']['f290lp']['footprint'] = []
        self.FileMap['NIRSPEC']['g395m']['f290lp']['file'] = []

        self.FileMap['NIRSPEC']['g395h'] = {}
        self.FileMap['NIRSPEC']['g395h']['f290lp'] = {}
        self.FileMap['NIRSPEC']['g395h']['f290lp']['footprint'] = []
        self.FileMap['NIRSPEC']['g395h']['f290lp']['file'] = []

#********************************************************************************
    def set_file_table(self,
                       input_models,
                       input_filenames):
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
                        self.FileMap['MIRI'][channel[k]][subchannel]['file'].append(input_model)
                        self.FileMap['MIRI'][channel[k]][subchannel]['footprint'].append(None)
            #________________________________________________________________________________
            #NIRSPEC instrument
            #________________________________________________________________________________
                elif instrument == 'NIRSPEC':
                    fwa = input_model.meta.instrument.filter.lower()
                    gwa = input_model.meta.instrument.grating.lower()

                    self.FileMap['NIRSPEC'][gwa][fwa]['file'].append(input_model)
                    self.FileMap['NIRSPEC'][gwa][fwa]['footprint'].append(None)
                    instrument = input_model.meta.instrument.name.lower()
                    mod = importlib.import_module('.' + instrument, 'jwst.assign_wcs')
                else:

                    log.info('Instrument not valid for cube')
        return instrument, detector

class ErrorNoAssignWCS(Exception):
    pass
