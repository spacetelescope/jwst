from ..stpipe import Step
from .. import datamodels
from stcal.dark_current import dark_sub


__all__ = ["DarkCurrentStep"]


class DarkCurrentStep(Step):
    """
    DarkCurrentStep: Performs dark current correction by subtracting
    dark current reference data from the input science data model.
    """

    class_alias = "dark_current"

    spec = """
        dark_output = output_file(default = None) # Dark model or averaged dark subtracted
    """

    reference_file_types = ['dark']

    def process(self, input):

        # Open the input data model
        with datamodels.RampModel(input) as input_model:

            # Get the name of the dark reference file to use
            self.dark_name = self.get_reference_file(input_model, 'dark')
            self.log.info('Using DARK reference file %s', self.dark_name)

            # Check for a valid reference file
            if self.dark_name == 'N/A':
                self.log.warning('No DARK reference file found')
                self.log.warning('Dark current step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.dark = 'SKIPPED'
                return result

            # Create name for the intermediate dark, if desired.
            dark_output = self.dark_output
            if dark_output is not None:
                dark_output = self.make_output_path(
                    basepath=dark_output,
                    suffix=False
                )

            # Open the dark ref file data model - based on Instrument
            instrument = input_model.meta.instrument.name
            if instrument == 'MIRI':
                dark_model = datamodels.DarkMIRIModel(self.dark_name)
            else:
                dark_model = datamodels.DarkModel(self.dark_name)

            # Do the dark correction
            result = dark_sub.do_correction(
                input_model, dark_model, dark_output
            )

            out_data, dark_data = result

            if dark_data is not None and dark_data.save:
                save_dark_data_as_dark_model(dark_data, dark_model, instrument)
            dark_model.close()

            out_ramp = dark_output_data_2_ramp_model(out_data, input_model)

        return out_ramp


def save_dark_data_as_dark_model(dark_data, dark_model, instrument):
    """
    Save dark data from the dark current step as the appropriate dark model.

    Parameters
    ----------
    dark_data: DarkData
        Dark data used in the dark current step.

    dark_model: DarkMIRIModel or DarkModel
        The input dark model from reference.

    instrument: str
        The instrument name.
    """
    if instrument == "MIRI":
        out_dark_model = datamodels.DarkMIRIModel(
            data=dark_data.data,
            dq=dark_data.groupdq,
            err=dark_data.err)
    else:
        out_dark_model = datamodels.DarkModel(
            data=dark_data.data,
            dq=dark_data.groupdq,
            err=dark_data.err)
    out_dark_model.update(dark_model)

    out_dark_model.meta.exposure.nframes = dark_data.exp_nframes
    out_dark_model.meta.exposure.ngroups = dark_data.exp_ngroups
    out_dark_model.meta.exposure.groupgap = dark_data.exp_groupgap
    out_dark_model.save(dark_data.output_name)
    out_dark_model.close()


def dark_output_data_2_ramp_model(out_data, input_model):
    """
    Convert computed output data from the dark step to a RampModel.

    Parameters
    ----------
    out_data: ScienceData
        Computed science data from the dark current step.

    input_model: RampModel
        The input ramp model from which to subtract the dark current.

    Return
    ------
    out_model: RampModel
        The output ramp model from the dark current step.
    """

    if out_data.cal_step == "SKIPPED":
        # If processing was skipped in the lower-level routines,
        # just return the unmodified input model
        input_model.meta.cal_step.dark_sub = "SKIPPED"
        return input_model
    else:
        out_model = input_model.copy()
        out_model.meta.cal_step.dark_sub = out_data.cal_step
        out_model.data = out_data.data
        out_model.groupdq = out_data.groupdq
        out_model.pixeldq = out_data.pixeldq
        out_model.err = out_data.err
        return out_model
