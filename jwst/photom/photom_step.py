#! /usr/bin/env python

from ..stpipe import Step, cmdline
from . import photom
from .. import datamodels


class PhotomStep(Step):
    """
    PhotomStep: Module for loading photometric conversion infomation from
        reference files and attaching or applying them to the input science
        data model
    """

    reference_file_types = ['photom', 'area']

    def process(self, input_file):

        try:
            dm = datamodels.open(input_file)
        except IOError:
            self.log.error('Input can not be opened as a Model.')

        # Report the detected type of input model
        model_type = dm.__class__.__name__
        self.log.debug("Input is {}".format(model_type))
        if model_type not in ('CubeModel', 'ImageModel',
                              'IFUImageModel', 'MultiSlitModel'):
            self.log.warning("Input is not one of the supported model types: "
                             "CubeModel, ImageModel IFUImageModel or MultiSlitModel.")

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
        phot = photom.DataSet(dm, phot_filename, area_filename)
        output_obj = phot.do_all()

        output_obj.meta.cal_step.photom = 'COMPLETE'

        return output_obj
