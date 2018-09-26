from ..stpipe import Step
from .. import datamodels
from . import bar_shadow


__all__ = ["BarShadowStep"]


class BarShadowStep(Step):
    """
    BarShadowStep: Inserts the bar shadow and wavelength arrays
    into the data.

    Bar shadow correction depends on the position of a pixel along the slit
    and the wavelength. It is only applied to uniform sources and only for
    NRS MSA data.

    """

    spec = """
    """

    reference_file_types = ['barshadow']

    def process(self, input):

        # Open the input data model
        with datamodels.open(input) as input_model:
            exp_type = input_model.meta.exposure.type
            if exp_type == 'NRS_MSASPEC':
                # Get the name of the bar shadow reference file to use
                self.barshadow_name = self.get_reference_file(input_model,
                                                              'barshadow')
                self.log.info('Using BARSHADOW reference file %s',
                              self.barshadow_name)

                # Check for a valid reference file
                if self.barshadow_name == 'N/A':
                    self.log.warning('No BARSHADOW reference file found')
                    self.log.warning('Bar shadow step will be skipped')
                    result = input_model.copy()
                    result.meta.cal_step.barshadow = 'SKIPPED'
                    return result

                # Open the barshadow ref file data model
                barshadow_model = datamodels.BarshadowModel(self.barshadow_name)

                # Do the bar shadow correction
                result = bar_shadow.do_correction(input_model, barshadow_model)

                barshadow_model.close()
                result.meta.cal_step.barshadow = 'COMPLETE'
            else:
                input_model.meta.cal_step.barshadow = 'SKIPPED'
                result = input_model
        return result
