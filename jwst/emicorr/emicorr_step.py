#! /usr/bin/env python

from stdatamodels.jwst import datamodels
from jwst.stpipe import Step
from jwst.emicorr import emicorr


__all__ = ["EmiCorrStep"]


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
        user_supplied_reffile = string(default=None)  # ASDF user-supplied reference file
        save_intermediate_results = boolean(default=False)  # If True and a reference file is created on the fly, save it to disk
        skip = boolean(default=True)  # Skip the step
    """  # noqa: E501

    reference_file_types = ["emicorr"]

    def process(self, step_input):
        """
        Apply EMI correction to input data.

        Parameters
        ----------
        step_input : str or DataModel
            Input ramp data.

        Returns
        -------
        RampModel
            EMI corrected output datamodel.
        """
        # Open the input data model as a RampModel
        with datamodels.RampModel(step_input) as input_model:
            # Catch the cases to skip
            instrument = input_model.meta.instrument.name
            if instrument != "MIRI":
                self.log.warning(f"EMI correction not implemented for instrument: {instrument}")
                input_model.meta.cal_step.emicorr = "SKIPPED"
                return input_model

            readpatt = input_model.meta.exposure.readpatt
            allowed_readpatts = ["FAST", "FASTR1", "SLOW", "SLOWR1"]
            if readpatt.upper() not in allowed_readpatts:
                self.log.warning(f"EMI correction not implemented for read pattern: {readpatt}")
                input_model.meta.cal_step.emicorr = "SKIPPED"
                return input_model

            # Work on a copy
            result = input_model.copy()

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
                self.log.info("Correcting with reference file created on-the-fly.")

            elif self.user_supplied_reffile is None:
                emicorr_ref_filename = self.get_reference_file(result, "emicorr")
                # Skip the step if no reference file is found
                if emicorr_ref_filename == "N/A":
                    self.log.warning("No reference file found.")
                    self.log.warning("EMICORR step will be skipped")
                    result.meta.cal_step.emicorr = "SKIPPED"
                    return result
                else:
                    self.log.info(f"Using CRDS reference file: {emicorr_ref_filename}")
                    emicorr_model = datamodels.EmiModel(emicorr_ref_filename)

            else:
                self.log.info(f"Using user-supplied reference file: {self.user_supplied_reffile}")
                emicorr_model = datamodels.EmiModel(self.user_supplied_reffile)

            # Do the correction
            save_onthefly_reffile = None
            if self.save_intermediate_results:
                if emicorr_ref_filename is None and self.user_supplied_reffile is None:
                    # get the same full path as input file to save the on-the-fly reference file
                    emicorr_ref_filename = self.make_output_path(suffix="emi_ref_waves")
                    save_onthefly_reffile = emicorr_ref_filename
            emicorr_output = emicorr.apply_emicorr(
                result, emicorr_model, save_onthefly_reffile=save_onthefly_reffile, **pars
            )
            if emicorr_output is None:
                self.log.warning("Step skipped")
                result.meta.cal_step.emicorr = "SKIPPED"
                return result
            else:
                result = emicorr_output
                result.meta.cal_step.emicorr = "COMPLETE"

            # Cleanup
            del emicorr_model

        return result
