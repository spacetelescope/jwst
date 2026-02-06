import logging
import warnings

from stdatamodels.exceptions import ValidationWarning
from stdatamodels.jwst import datamodels

from jwst.emicorr import emicorr
from jwst.stpipe import Step

__all__ = ["EmiCorrStep"]

log = logging.getLogger(__name__)


class EmiCorrStep(Step):
    """Correct MIRI ramp data for EMI noise."""

    class_alias = "emicorr"

    spec = """
        algorithm = option('sequential', 'joint', default='joint')  # EMI fitting algorithm
        nints_to_phase = integer(default=None)  # Number of integrations to phase
        nbins = integer(default=None)  # Number of bins in one phased wave
        scale_reference = boolean(default=True)  # If True, the reference wavelength will be scaled to the data's phase amplitude
        onthefly_corr_freq = float_list(default=None)  # Frequencies to use for correction
        use_n_cycles = integer(default=3)  # Use N cycles to calculate the phase, to use all integrations set to None
        fit_ints_separately = boolean(default=False)  # If True and algorithm is 'joint', each integration is separately fit.
        user_supplied_reffile = string(default=None)  # Deprecated; use 'override_emicorr' instead
        save_intermediate_results = boolean(default=False)  # If True and a reference file is created on the fly, save it to disk
        skip = boolean(default=True)  # Skip the step
    """  # noqa: E501

    reference_file_types = ["emicorr"]

    def process(self, step_input):
        """
        Apply EMI correction to input data.

        Parameters
        ----------
        step_input : str or `~stdatamodels.jwst.datamodels.RampModel`
            Input ramp filename or datamodel.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.RampModel`
            EMI corrected output datamodel.
        """
        # Open the input data model as a RampModel, catching a specific
        # expected warning for uncal files
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "error",
                    message=r"(?s:.*)Array datatype .* not compatible",
                    category=ValidationWarning,
                )
                result = self.prepare_output(step_input, open_as_type=datamodels.RampModel)
        except ValidationWarning as err:
            log.error(err)

            # Inform the user and raise a clearer error message
            msg = (
                "Input data model does not have float-type data. "
                "The file should be opened as a RampModel before calling the step."
            )
            log.error(msg)
            raise TypeError(msg) from None

        # Catch the cases to skip
        instrument = result.meta.instrument.name
        if instrument != "MIRI":
            log.warning(f"EMI correction not implemented for instrument: {instrument}")
            result.meta.cal_step.emicorr = "SKIPPED"
            return result

        readpatt = result.meta.exposure.readpatt
        allowed_readpatts = ["FAST", "FASTR1", "SLOW", "SLOWR1"]
        if readpatt.upper() not in allowed_readpatts:
            log.warning(f"EMI correction not implemented for read pattern: {readpatt}")
            result.meta.cal_step.emicorr = "SKIPPED"
            return result

        # Setup parameters
        pars = {
            "algorithm": self.algorithm,
            "nints_to_phase": self.nints_to_phase,
            "nbins": self.nbins,
            "scale_reference": self.scale_reference,
            "onthefly_corr_freq": self.onthefly_corr_freq,
            "use_n_cycles": self.use_n_cycles,
            "fit_ints_separately": self.fit_ints_separately,
        }

        # Get the reference file
        save_onthefly_reffile, emicorr_ref_filename, emicorr_model = None, None, None
        if self.onthefly_corr_freq is not None:
            emicorr_ref_filename = None
            log.info("Correcting with reference file created on-the-fly.")

        else:
            # Check for a deprecated parameter
            if self.user_supplied_reffile is not None:
                msg = (
                    "The 'user_supplied_reffile' parameter is deprecated and will "
                    "be removed in a future release. Use 'override_emicorr' instead."
                )
                warnings.warn(msg, DeprecationWarning, stacklevel=2)
                log.warning(msg)
                emicorr_ref_filename = self.user_supplied_reffile
            else:
                emicorr_ref_filename = self.get_reference_file(result, "emicorr")

            # Skip the step if no reference file is found
            if emicorr_ref_filename == "N/A":
                log.warning("No reference file found.")
                log.warning("EMICORR step will be skipped")
                result.meta.cal_step.emicorr = "SKIPPED"
                return result
            else:
                log.info(f"Using EMICORR reference file: {emicorr_ref_filename}")
                emicorr_model = datamodels.EmiModel(emicorr_ref_filename)

        # Do the correction
        save_onthefly_reffile = None
        if self.save_intermediate_results:
            if emicorr_ref_filename is None:
                # get the same full path as input file to save the on-the-fly reference file
                emicorr_ref_filename = self.make_output_path(suffix="emi_ref_waves")
                save_onthefly_reffile = emicorr_ref_filename
        emicorr_output = emicorr.apply_emicorr(
            result, emicorr_model, save_onthefly_reffile=save_onthefly_reffile, **pars
        )
        if emicorr_output is None:
            log.warning("Step skipped")
            result.meta.cal_step.emicorr = "SKIPPED"
            return result
        else:
            result = emicorr_output
            result.meta.cal_step.emicorr = "COMPLETE"

        # Cleanup
        del emicorr_model

        return result
