2#! /usr/bin/env python

from jwst.stpipe import Step, cmdline
from .. import datamodels
from . import flat_field


class FlatFieldStep(Step):
    """
    FlatFieldStep: Flat-field a science image using a flatfield reference image.
    """

    spec = """
        # Suffix for optional output file for interpolated flat fields.
        flat_suffix = string(default=None)
    """

    reference_file_types = ['flat']

    def process(self, input):

        input_model = datamodels.open(input)

        # Figure out what kind of input data model is in use
        if isinstance(input_model, datamodels.CubeModel): # multi-integration dataset
            self.log.debug('Input is a CubeModel')
        elif isinstance(input_model, datamodels.ImageModel):
            self.log.debug('Input is an ImageModel')
        elif isinstance(input_model, datamodels.DataModel):
            self.log.debug('Input is a MultiSlitModel')
            input_model.close()
            input_model = datamodels.MultiSlitModel(input)

        # Retrieve the reference file name
        self.flat_filename = self.get_reference_file(input_model, 'flat')
        self.log.info('Using FLAT reference file: %s', self.flat_filename)

        # Check for a valid reference file
        if self.flat_filename == 'N/A':
            self.log.warning('No FLAT reference file found')
            self.log.warning('Flat-field step will be skipped')
            result = input_model.copy()
            result.meta.cal_step.flat_field = 'SKIPPED'
            input_model.close()
            return result

        # Find out what model to use for the flat field reference file.
        # If datamodels.open() thinks it's a CubeModel, open it as a
        # CubeFlatModel; otherwise, open it as a MultiSlitModel.
        flat_model = datamodels.open(self.flat_filename)
        if isinstance(flat_model, datamodels.CubeModel):
            # MOS/MSA_mode
            self.log.debug('Flatfield is a CubeFlatModel')
            flat_model.close()
            flat_model = datamodels.CubeFlatModel(self.flat_filename)
        else:
            self.log.debug('Flatfield is an ImageModel or MultiSlitModel')
            flat_model.close()
            flat_model = datamodels.MultiSlitModel(self.flat_filename)

        # Do the flat-field correction
        (output_model, interpolated_flats) = \
                flat_field.do_correction(input_model, flat_model,
                                         self.flat_suffix)

        # Close the inputs
        input_model.close()
        flat_model.close()

        if interpolated_flats is not None:
            self.log.info("Writing interpolated flat fields.")
            self.save_model(interpolated_flats, self.flat_suffix)
            interpolated_flats.close()

        return output_model

if __name__ == '__main__':
    cmdline.step_script(flat_field_step)
