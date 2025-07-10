from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from stcal.dark_current import dark_sub
import numpy as np


__all__ = ["DarkCurrentStep"]


class DarkCurrentStep(Step):
    """Subtract dark current from the input science data."""

    class_alias = "dark_current"

    spec = """
        dark_output = output_file(default = None) # Dark model or averaged dark subtracted
        average_dark_current = float(default=None) # The average dark current for this detector in units of e-/sec.
    """  # noqa: E501

    reference_file_types = ["dark"]

    def process(self, step_input):
        """
        Perform the dark current subtraction step.

        Returns
        -------
        out_ramp : `~stdatamodels.jwst.datamodels.RampModel`
            The science model with dark current subtracted form the data array.
        """
        # Open the input data model
        with datamodels.RampModel(step_input) as input_model:
            # Get the name of the dark reference file to use
            self.dark_name = self.get_reference_file(input_model, "dark")
            self.log.info("Using DARK reference file %s", self.dark_name)

            # Check for a valid reference file
            if self.dark_name == "N/A":
                self.log.warning("No DARK reference file found")
                self.log.warning("Dark current step will be skipped")
                input_model.meta.cal_step.dark = "SKIPPED"
                return input_model

            # Work on a copy
            result = input_model.copy()

            # Create name for the intermediate dark, if desired.
            dark_output = self.dark_output
            if dark_output is not None:
                dark_output = self.make_output_path(basepath=dark_output, suffix=False)

            # Open the dark ref file data model - based on Instrument
            instrument = result.meta.instrument.name
            if instrument == "MIRI":
                dark_model = datamodels.DarkMIRIModel(self.dark_name)
            elif instrument == "NIRSPEC":
                dark_model = datamodels.DarkNirspecModel(self.dark_name)
            else:
                dark_model = datamodels.DarkModel(self.dark_name)

            # Store user-defined average_dark_current in model, if provided
            # A user-defined value will take precedence over any value present
            # in dark reference file
            self.set_average_dark_current(result, dark_model)

            # Do the dark correction
            correction = dark_sub.do_correction(result, dark_model, dark_output)

            out_data, dark_data = correction

            if dark_data is not None and dark_data.save:
                save_dark_data_as_dark_model(dark_data, dark_model, instrument)

            out_ramp = dark_output_data_2_ramp_model(out_data, result)

            # Cleanup
            del dark_model
            del result

        return out_ramp

    def set_average_dark_current(self, input_model, dark_model):
        """
        Assign average dark current.

        Take the three possible locations specifying
        the average dark current and assign them to the
        input model, in priority order:
        1) Any value provided to the step parameter, either from
        the user or a parameter reference file
        2) The 2-D array stored in dark_model.average_dark_current
        3) The scalar value stored in dark_model.meta.exposure.average_dark_current

        Parameters
        ----------
        input_model : `~stdatamodels.jwst.datamodels.RampModel`
            The input datamodel containing the 4-D ramp array.
        dark_model : Union[stdatamodels.jwst.datamodels.DarkModel,
                           stdatamodels.jwst.datamodels.DarkMIRIModel]
            The dark reference file datamodel.
        """
        if self.average_dark_current is not None:
            input_model.average_dark_current[:, :] = self.average_dark_current
            self.log.info(
                "Using Poisson noise from average dark current %s e-/sec", self.average_dark_current
            )
        else:
            # First prioritize a 2D average_dark_current, if defined in dark.
            # If not present, apply scalar to datamodel array, if scalar is present.
            if np.sum(dark_model.average_dark_current) == 0.0:
                input_model.average_dark_current[:, :] = (
                    dark_model.meta.exposure.average_dark_current
                )
            elif np.shape(input_model.average_dark_current) != np.shape(
                dark_model.average_dark_current
            ):
                self.log.warning("DarkModel average_dark_current does not match shape of data.")
                self.log.warning("Dark current from reference file cannot be applied.")
            else:
                input_model.average_dark_current = dark_model.average_dark_current


def save_dark_data_as_dark_model(dark_data, dark_model, instrument):
    """
    Save dark data from the dark current step as the appropriate dark model.

    Parameters
    ----------
    dark_data : ndarray
        Dark data used in the dark current step.

    dark_model : `~stdatamodels.jwst.datamodels.DarkMIRIModel` or
                 `~stdatamodels.jwst.datamodels.DarkModel`
        The input dark model from reference.

    instrument : str
        The instrument name.
    """
    if instrument == "MIRI":
        out_dark_model = datamodels.DarkMIRIModel(data=dark_data.data, dq=dark_data.groupdq)
    else:
        out_dark_model = datamodels.DarkModel(data=dark_data.data, dq=dark_data.groupdq)
    out_dark_model.update(dark_model)

    out_dark_model.meta.exposure.nframes = dark_data.exp_nframes
    out_dark_model.meta.exposure.ngroups = dark_data.exp_ngroups
    out_dark_model.meta.exposure.groupgap = dark_data.exp_groupgap
    out_dark_model.save(dark_data.output_name)
    out_dark_model.close()


def dark_output_data_2_ramp_model(out_data, out_model):
    """
    Convert computed output data from the dark step to a RampModel.

    Parameters
    ----------
    out_data : `~stdatamodels.jwst.datamodels.DataModel`
        Computed science data from the dark current step.

    out_model : `~stdatamodels.jwst.datamodels.RampModel`
        The input ramp model from which to subtract the dark current.

    Returns
    -------
    out_model : `~stdatamodels.jwst.datamodels.RampModel`
        The output ramp model from the dark current step.
    """
    if out_data.cal_step == "SKIPPED":
        # If processing was skipped in the lower-level routines,
        # just return the unmodified input model
        out_model.meta.cal_step.dark_sub = "SKIPPED"
        return out_model
    else:
        out_model.meta.cal_step.dark_sub = out_data.cal_step
        out_model.data = out_data.data
        out_model.groupdq = out_data.groupdq
        out_model.pixeldq = out_data.pixeldq
        return out_model
