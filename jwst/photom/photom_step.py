#! /usr/bin/env python

from ..stpipe import Step
from . import photom
from .. import datamodels

__all__ = ["PhotomStep"]


class PhotomStep(Step):
    """
    PhotomStep: Module for loading photometric conversion infomation from
        reference files and attaching or applying them to the input science
        data model
    """

    reference_file_types = ['photom', 'area']

    def process(self, input):

        try:
            input_model = datamodels.open(input)
        except IOError:
            self.log.error('Input can not be opened as a Model.')

        # Report the detected type of input model
        model_type = input_model.__class__.__name__
        self.log.debug("Input is {}".format(model_type))
        if model_type not in ('CubeModel', 'ImageModel', 'SlitModel',
                              'IFUImageModel', 'MultiSlitModel'):
            self.log.warning("Input is not one of the supported model types: "
                             "CubeModel, ImageModel, IFUImageModel or "
                             "MultiSlitModel.")

        # Get the reference file names
        phot_filename = self.get_reference_file(input_model, 'photom')
        self.log.info('Using photom reference file: %s', phot_filename)
        area_filename = self.get_reference_file(input_model, 'area')
        self.log.info('Using area reference file: %s', area_filename)

        # Check for a valid photom reference file
        if phot_filename == 'N/A':
            self.log.warning('No PHOTOM reference file found')
            self.log.warning('Photom step will be skipped')
            result = input_model.copy()
            result.meta.cal_step.photom = 'SKIPPED'
            return result

        # Do the correction
        phot = photom.DataSet(input_model)
        result = phot.apply_photom(phot_filename, area_filename)

        result.meta.cal_step.photom = 'COMPLETE'

        return result
