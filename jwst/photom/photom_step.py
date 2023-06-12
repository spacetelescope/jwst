#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import photom

__all__ = ["PhotomStep"]


class PhotomStep(Step):
    """
    PhotomStep: Module for loading photometric conversion information from
        reference files and attaching or applying them to the input science
        data model
    """

    class_alias = "photom"

    spec = """
        inverse = boolean(default=False)    # Invert the operation
        source_type = string(default=None)  # Process as specified source type.
        mrs_time_correction = boolean(default=True) # Apply the MIRI MRS time dependent correction
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
                              'IFUImageModel', 'MultiSlitModel',
                              'MultiSpecModel'):
            self.log.warning("Input is not one of the supported model types: "
                             "CubeModel, ImageModel, IFUImageModel, "
                             "MultiSlitModel, or MultiSpecModel.")

        # Setup reference files and whether previous correction information
        # should be used.
        if self.use_correction_pars and self.correction_pars:
            self.log.info('Using previously specified correction parameters.')
            correction_pars = self.correction_pars
            phot_filename = correction_pars['refs']['photom']
            area_filename = correction_pars['refs']['area']
        else:
            correction_pars = None
            phot_filename = self.get_reference_file(input_model, 'photom')
            area_filename = self.get_reference_file(input_model, 'area')

        self.log.info('Using photom reference file: %s', phot_filename)
        self.log.info('Using area reference file: %s', area_filename)

        # Check for a valid photom reference file
        if phot_filename == 'N/A':
            self.log.warning('No PHOTOM reference file found')
            self.log.warning('Photom step will be skipped')
            result = input_model.copy()
            result.meta.cal_step.photom = 'SKIPPED'
            return result

        # Do the correction
        phot = photom.DataSet(input_model, self.inverse, self.source_type,
                              self.mrs_time_correction, correction_pars)
        result = phot.apply_photom(phot_filename, area_filename)
        result.meta.cal_step.photom = 'COMPLETE'

        self.correction_pars = phot.correction_pars
        self.correction_pars['refs'] = {'photom': phot_filename, 'area': area_filename}

        return result
