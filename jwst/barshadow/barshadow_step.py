from stdatamodels.jwst import datamodels

from ..stpipe import Step
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

    class_alias = "barshadow"

    spec = """
        inverse = boolean(default=False)    # Invert the operation
        source_type = string(default=None)  # Process as specified source type.
    """

    reference_file_types = ['barshadow']

    def process(self, input):
        """Perform the barshadow correction step

        Parameters
        ----------
        input : JWST datamodel
            input JWST datamodel object

        Returns
        -------
        result : jwst datamodel
            JWST datamodel object with barshadow extension(s) added
        """

        # Open the input data model
        with datamodels.open(input) as input_model:
            exp_type = input_model.meta.exposure.type
            if exp_type == 'NRS_MSASPEC':

                if self.use_correction_pars:
                    correction_pars = self.correction_pars
                    barshadow_model = None

                else:
                    correction_pars = None

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
                result, self.correction_pars = bar_shadow.do_correction(
                    input_model, barshadow_model,
                    inverse=self.inverse, source_type=self.source_type,
                    correction_pars=correction_pars
                )

                if barshadow_model:
                    barshadow_model.close()
                result.meta.cal_step.barshadow = 'COMPLETE'
            else:
                input_model.meta.cal_step.barshadow = 'SKIPPED'
                result = input_model
        return result
