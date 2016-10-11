from ..stpipe import Step
from .. import datamodels
from . import path_loss


class PathLossStep(Step):
    """
    PathLossStep: Performs path loss correction by interpolating
    path loss reference data over wavlength for each region.
    """

    spec = """
    """

    reference_file_types = ['pathloss']

    def process(self, input):

        # Open the input data model
        with datamodels.open(input) as input_model:

            # Get the name of the path loss reference file to use
            self.pathloss_name = self.get_reference_file(input_model,
                                                         'pathloss')
            self.log.info('Using PATHLOSS reference file %s',
                          self.pathloss_name)

            # Check for a valid reference file
            if self.pathloss_name == 'N/A':
                self.log.warning('No PATHLOSS reference file found')
                self.log.warning('Path loss step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.pathloss = 'SKIPPED'
                return result

            instrument = input_model.meta.instrument.name
            # Open the pathloss ref file data model
            pathloss_model = datamodels.PathlossModel(self.pathloss_name)

            # Do the pathloss correction     
            result = path_loss.do_correction(input_model, pathloss_model)

            pathloss_model.close()

        return result
