""" Dictionary holding defaults for cube_build
"""
from stdatamodels.jwst import datamodels
from jwst.datamodels import ModelContainer
import logging
from os import path
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class FileTable():
    """ Dictionary contains defaults for MIRI and NIRSPEC data
    """

    def __init__(self, in_memory):

        self.in_memory = in_memory
        self.FileMap = {}
        self.FileMap['filename'] = []
        self.FileMap['MIRI'] = {}
        
        self.FileMap['MIRI']['1'] = {}
        self.FileMap['MIRI']['1']['short'] = []
        self.FileMap['MIRI']['1']['medium'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['1']['long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['1']['short-medium'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['1']['short-long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['1']['medium-short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['1']['medium-long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['1']['long-short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['1']['long-medium'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['MIRI']['2'] = {}
        self.FileMap['MIRI']['2']['short'] = []
        self.FileMap['MIRI']['2']['medium'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['2']['long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['2']['short-medium'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['2']['short-long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['2']['medium-short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['2']['medium-long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['2']['long-short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['2']['long-medium'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['MIRI']['3'] = {}
        self.FileMap['MIRI']['3']['short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['3']['medium'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['3']['long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['3']['short-medium'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['3']['short-long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['3']['medium-short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['3']['medium-long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['3']['long-short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['3']['long-medium'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['MIRI']['4'] = {}
        self.FileMap['MIRI']['4']['short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['4']['medium'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['4']['long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['4']['short-medium'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['4']['short-long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['4']['medium-short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['4']['medium-long'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['4']['long-short'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['MIRI']['4']['long-medium'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['NIRSPEC'] = {}
        self.FileMap['NIRSPEC']['prism'] = {}
        self.FileMap['NIRSPEC']['prism']['clear'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['NIRSPEC']['prism']['opaque'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['NIRSPEC']['g140m'] = {}
        self.FileMap['NIRSPEC']['g140m']['f070lp'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['NIRSPEC']['g140m']['f100lp'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['NIRSPEC']['g140m']['opaque'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['NIRSPEC']['g140h'] = {}
        self.FileMap['NIRSPEC']['g140h']['f070lp'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['NIRSPEC']['g140h']['f100lp'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['NIRSPEC']['g140h']['opaque'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['NIRSPEC']['g235m'] = {}
        self.FileMap['NIRSPEC']['g235m']['f170lp'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['NIRSPEC']['g235m']['opaque'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['NIRSPEC']['g235h'] = {}
        self.FileMap['NIRSPEC']['g235h']['f170lp'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['NIRSPEC']['g235h']['opaque'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['NIRSPEC']['g395m'] = {}
        self.FileMap['NIRSPEC']['g395m']['f290lp'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['NIRSPEC']['g395m']['opaque'] = ModelContainer(save_open=self.in_memory)

        self.FileMap['NIRSPEC']['g395h'] = {}
        self.FileMap['NIRSPEC']['g395h']['f290lp'] = ModelContainer(save_open=self.in_memory)
        self.FileMap['NIRSPEC']['g395h']['opaque'] = ModelContainer(save_open=self.in_memory)

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

        print('in file table',input_models._models)
        if not self.in_memory:
            catfile_dir = path.dirname(input_models._models[0])
            print('catfile_dir',catfile_dir)
        for model in input_models:
            instrument = model.meta.instrument.name.upper()
            assign_wcs = model.meta.cal_step.assign_wcs
            
            filename  = model.meta.filename
            print(filename)
            self.FileMap['filename'].append(filename)

            if not isinstance(model, datamodels.IFUImageModel):
                        raise NotIFUImageModel(
                           f"Input data is not a IFUImageModel, instead it is {model}")    
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
                    if self.in_memory:
                        self.FileMap['MIRI'][channel[k]][subchannel].append(model.copy())
                    else:
                        self.FileMap['MIRI'][channel[k]][subchannel].append(catfile_dir + '/' + filename)   
            # _____________________________________________________________________
            # NIRSPEC instrument
            elif instrument == 'NIRSPEC':
                fwa = model.meta.instrument.filter.lower()
                gwa = model.meta.instrument.grating.lower()
                self.FileMap['NIRSPEC'][gwa][fwa].append(model.copy())
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
