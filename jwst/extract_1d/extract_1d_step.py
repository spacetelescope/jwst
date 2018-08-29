#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import extract


__all__ = ["Extract1dStep"]


class Extract1dStep(Step):
    """
    Extract1dStep: Extract a 1-d spectrum from 2-d data
    """

    spec = """
    # Boxcar smoothing width for background regions.
    smoothing_length = integer(default=None)
    # Order of polynomial fit to one column (or row if the dispersion
    # direction is vertical) of background regions.
    bkg_order = integer(default=None, min=0)
    # Log a progress message when processing multi-integration data.
    log_increment = integer(default=50)
    """

    reference_file_types = ['extract1d']

    def process(self, input):

        # Open the input and figure out what type of model it is
        input_model = datamodels.open(input)

        if isinstance(input_model, datamodels.CubeModel):
            # It's a 3-D multi-integration model
            self.log.debug('Input is a CubeModel for a multiple integ. file')
        elif isinstance(input_model, datamodels.ImageModel):
            # It's a single 2-D image
            self.log.debug('Input is an ImageModel')
        elif isinstance(input_model, datamodels.ModelContainer):
            self.log.debug('Input is a ModelContainer')
        elif isinstance(input_model, datamodels.MultiSlitModel):
            self.log.debug('Input is a MultiSlitModel')
        elif isinstance(input_model, datamodels.MultiProductModel):
            self.log.debug('Input is a MultiProductModel')
        elif isinstance(input_model, datamodels.IFUCubeModel):
            self.log.debug('Input is an IFUCubeModel')
        elif isinstance(input_model, datamodels.DrizProductModel):
            # Resampled 2-D data
            self.log.debug('Input is a DrizProductModel')
        elif isinstance(input_model, datamodels.SlitModel):
            # NRS_BRIGHTOBJ mode
            self.log.debug('Input is a SlitModel')
        else:
            self.log.error('Input is a %s,', str(type(input_model)))
            self.log.error('which was not expected for extract_1d.')
            self.log.error('extract_1d will be skipped.')
            input_model.meta.cal_step.extract_1d = 'SKIPPED'
            return input_model

        # Do the extraction
        if isinstance(input_model, datamodels.ModelContainer):
            if len(input_model) > 1:
                self.log.debug("Input contains %d items", len(input_model))
                result = datamodels.ModelContainer()
                for model in input_model:
                    if model.meta.exposure.type in extract.WFSS_EXPTYPES:
                        self.ref_file = 'N/A'
                        self.log.info('No EXTRACT1D reference file '
                                      'will be used')
                    else:
                        # Get the reference file name
                        self.ref_file = self.get_reference_file(
                                        model, 'extract1d')
                        self.log.info('Using EXTRACT1D reference file %s',
                                      self.ref_file)
                    temp = extract.do_extract1d(model, self.ref_file,
                                                self.smoothing_length,
                                                self.bkg_order,
                                                self.log_increment)
                    # Set the step flag to complete in each MultiSpecModel
                    temp.meta.cal_step.extract_1d = 'COMPLETE'
                    result.append(temp)
                    del temp
            elif len(input_model) == 1:
                if input_model[0].meta.exposure.type in extract.WFSS_EXPTYPES:
                    self.ref_file = 'N/A'
                    self.log.info('No EXTRACT1D reference file will be used')
                else:
                    # Get the reference file name for the one model in input
                    self.ref_file = self.get_reference_file(input_model[0],
                                                            'extract1d')
                    self.log.info('Using EXTRACT1D reference file %s',
                                  self.ref_file)
                result = extract.do_extract1d(input_model[0], self.ref_file,
                                              self.smoothing_length,
                                              self.bkg_order,
                                              self.log_increment)
                # Set the step flag to complete
                result.meta.cal_step.extract_1d = 'COMPLETE'
            else:
                self.log.error('Input model is empty;')
                self.log.error('extract_1d will be skipped.')
                return input_model
        else:
            # Get the reference file name
            if input_model.meta.exposure.type in extract.WFSS_EXPTYPES:
                self.ref_file = 'N/A'
                self.log.info('No EXTRACT1D reference file will be used')
            else:
                self.ref_file = self.get_reference_file(input_model,
                                                        'extract1d')
                self.log.info('Using EXTRACT1D reference file %s',
                              self.ref_file)
            result = extract.do_extract1d(input_model, self.ref_file,
                                          self.smoothing_length,
                                          self.bkg_order,
                                          self.log_increment)
            # Set the step flag to complete
            result.meta.cal_step.extract_1d = 'COMPLETE'

        input_model.close()

        return result


def extract_1d_correction(input):
    a = Extract1dStep()
    result = a.process(input)
    return result
