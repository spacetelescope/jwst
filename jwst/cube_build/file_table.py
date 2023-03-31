""" Dictionary holding defaults for cube_build
"""
from stdatamodels.jwst import datamodels
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class FileTable():
    """ Dictionary contains defaults for MIRI and NIRSPEC data
    """

    def __init__(self):

        self.FileMap = {}
        self.FileMap['filename'] = []
        self.FileMap['MIRI'] = {}
        self.FileMap['MIRI']['1'] = {}
        self.FileMap['MIRI']['1']['short'] = []
        self.FileMap['MIRI']['1']['medium'] = []
        self.FileMap['MIRI']['1']['long'] = []
        self.FileMap['MIRI']['1']['short-medium'] = []
        self.FileMap['MIRI']['1']['short-long'] = []
        self.FileMap['MIRI']['1']['medium-short'] = []
        self.FileMap['MIRI']['1']['medium-long'] = []
        self.FileMap['MIRI']['1']['long-short'] = []
        self.FileMap['MIRI']['1']['long-medium'] = []

        self.FileMap['MIRI']['2'] = {}
        self.FileMap['MIRI']['2']['short'] = []
        self.FileMap['MIRI']['2']['medium'] = []
        self.FileMap['MIRI']['2']['long'] = []
        self.FileMap['MIRI']['2']['short-medium'] = []
        self.FileMap['MIRI']['2']['short-long'] = []
        self.FileMap['MIRI']['2']['medium-short'] = []
        self.FileMap['MIRI']['2']['medium-long'] = []
        self.FileMap['MIRI']['2']['long-short'] = []
        self.FileMap['MIRI']['2']['long-medium'] = []

        self.FileMap['MIRI']['3'] = {}
        self.FileMap['MIRI']['3']['short'] = []
        self.FileMap['MIRI']['3']['medium'] = []
        self.FileMap['MIRI']['3']['long'] = []
        self.FileMap['MIRI']['3']['short-medium'] = []
        self.FileMap['MIRI']['3']['short-long'] = []
        self.FileMap['MIRI']['3']['medium-short'] = []
        self.FileMap['MIRI']['3']['medium-long'] = []
        self.FileMap['MIRI']['3']['long-short'] = []
        self.FileMap['MIRI']['3']['long-medium'] = []

        self.FileMap['MIRI']['4'] = {}
        self.FileMap['MIRI']['4']['short'] = []
        self.FileMap['MIRI']['4']['medium'] = []
        self.FileMap['MIRI']['4']['long'] = []
        self.FileMap['MIRI']['4']['short-medium'] = []
        self.FileMap['MIRI']['4']['short-long'] = []
        self.FileMap['MIRI']['4']['medium-short'] = []
        self.FileMap['MIRI']['4']['medium-long'] = []
        self.FileMap['MIRI']['4']['long-short'] = []
        self.FileMap['MIRI']['4']['long-medium'] = []

        self.FileMap['NIRSPEC'] = {}
        self.FileMap['NIRSPEC']['prism'] = {}
        self.FileMap['NIRSPEC']['prism']['clear'] = []
        self.FileMap['NIRSPEC']['prism']['opaque'] = []

        self.FileMap['NIRSPEC']['g140m'] = {}
        self.FileMap['NIRSPEC']['g140m']['f070lp'] = []
        self.FileMap['NIRSPEC']['g140m']['f100lp'] = []
        self.FileMap['NIRSPEC']['g140m']['opaque'] = []

        self.FileMap['NIRSPEC']['g140h'] = {}
        self.FileMap['NIRSPEC']['g140h']['f070lp'] = []
        self.FileMap['NIRSPEC']['g140h']['f100lp'] = []
        self.FileMap['NIRSPEC']['g140h']['opaque'] = []

        self.FileMap['NIRSPEC']['g235m'] = {}
        self.FileMap['NIRSPEC']['g235m']['f170lp'] = []
        self.FileMap['NIRSPEC']['g235m']['opaque'] = []

        self.FileMap['NIRSPEC']['g235h'] = {}
        self.FileMap['NIRSPEC']['g235h']['f170lp'] = []
        self.FileMap['NIRSPEC']['g235h']['opaque'] = []

        self.FileMap['NIRSPEC']['g395m'] = {}
        self.FileMap['NIRSPEC']['g395m']['f290lp'] = []
        self.FileMap['NIRSPEC']['g395m']['opaque'] = []

        self.FileMap['NIRSPEC']['g395h'] = {}
        self.FileMap['NIRSPEC']['g395h']['f290lp'] = []
        self.FileMap['NIRSPEC']['g395h']['opaque'] = []

# ********************************************************************************

    def set_file_table(self,
                       input_models):
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
        instrument name
        """
# ________________________________________________________________________________
# Loop over input list of files and assign fill in the MasterTable with filename
# for the correct (channel-subchannel) or (grating-subchannel)

        for model in input_models:
            instrument = model.meta.instrument.name.upper()
            assign_wcs = model.meta.cal_step.assign_wcs
            filename = model.meta.filename
            self.FileMap['filename'].append(filename)

            if not isinstance(model, datamodels.IFUImageModel):
                raise NotIFUImageModel(f"Input data is not a IFUImageModel, instead it is {model}")
            if assign_wcs != 'COMPLETE':
                raise ErrorNoAssignWCS("Assign WCS has not been run on file %s",
                                       model.meta.filename)
            # _____________________________________________________________________
            # MIRI instrument
            if instrument == 'MIRI':
                channel = model.meta.instrument.channel
                subchannel = model.meta.instrument.band.lower()
                clenf = len(channel)
                for k in range(clenf):
                    self.FileMap['MIRI'][channel[k]][subchannel].append(model)
            # _____________________________________________________________________
            # NIRSPEC instrument
            elif instrument == 'NIRSPEC':
                fwa = model.meta.instrument.filter.lower()
                gwa = model.meta.instrument.grating.lower()
                self.FileMap['NIRSPEC'][gwa][fwa].append(model)
            else:
                log.info('Instrument not valid for cube')
                pass

            model.close()
            del model

        return instrument


class ErrorNoAssignWCS(Exception):
    pass


class NotIFUImageModel(Exception):
    """ Raise Exception if data is not of type IFUImageModel
    """
    pass
