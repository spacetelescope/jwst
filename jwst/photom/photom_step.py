#! /usr/bin/env python

from ..stpipe import Step, cmdline
from . import photom
from .. import datamodels


class PhotomStep(Step):
    """
    PhotomStep: Module for extraction photom conversion factor(s)
        and writing them to input header
    """

    reference_file_types = ['photom', 'area']

    def process(self, input_file):

        try:
            dm = datamodels.open(input_file)
        except IOError:
            self.log.error('Input can not be opened as a Model.')

        # Open input as correct type
        if isinstance(dm, datamodels.CubeModel): # _integ.fits product: 3D array
            self.log.debug('Input is a CubeModel for a multiple integ file.')
        elif isinstance(dm, datamodels.ImageModel):  # standard product: 2D array
            self.log.debug('Input is an ImageModel.')
        elif isinstance(dm, datamodels.DataModel): # multi 2D arrays
            self.log.debug('Input is a MultiSlitModel.')
            dm.close()
            dm = datamodels.MultiSlitModel(input_file)
        else:
            self.log.warning('Input is not a CubeModel, ImageModel or MultiSlitModel.')

        # Get the reference file names
        phot_filename = self.get_reference_file(dm, 'photom')
        self.log.info('Using photom reference file: %s', phot_filename)
        area_filename = self.get_reference_file(dm, 'area')
        self.log.info('Using area reference file: %s', area_filename)

        # Check for a valid photom reference file
        if phot_filename == 'N/A':
            self.log.warning('No PHOTOM reference file found')
            self.log.warning('Photom step will be skipped')
            result = dm.copy()
            result.meta.cal_step.photom = 'SKIPPED'
            return result

        # Do the correction
        ff_a = photom.DataSet(dm, phot_filename, area_filename)
        output_obj = ff_a.do_all()

        output_obj.meta.cal_step.photom = 'COMPLETE'

        return output_obj

if __name__ == '__main__':
    cmdline.step_script(photom_step)
