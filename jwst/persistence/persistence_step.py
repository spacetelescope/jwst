#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from . import persistence

__all__ = ["PersistenceStep"]


class PersistenceStep(Step):
    """PersistenceStep: Correct a science image for persistence."""

    class_alias = "persistence"

    spec = """
        input_trapsfilled = string(default="") # Name of the most recent trapsfilled file for the current detector
        flag_pers_cutoff = float(default=40.) # Pixels with persistence correction >= this value in DN will be flagged in the DQ
        save_persistence = boolean(default=False) # Save subtracted persistence to an output file with suffix '_output_pers'
        save_trapsfilled = boolean(default=True) # Save updated trapsfilled file with suffix '_trapsfilled'
        modify_input = boolean(default=False)
    """  # noqa: E501

    reference_file_types = ["trapdensity", "trappars", "persat"]

    def process(self, step_input):
        """
        Execute the persistence correction step.

        Parameters
        ----------
        step_input : DataModel or str
            Input datamodel or file to be corrected

        Returns
        -------
        output_model : DataModel
            The persistence corrected datamodel
        """
        if self.input_trapsfilled is not None:
            if (self.input_trapsfilled == "None") or (len(self.input_trapsfilled) == 0):
                self.input_trapsfilled = None

        with datamodels.RampModel(step_input) as input_model:
            self.trap_density_filename = self.get_reference_file(input_model, "trapdensity")
            self.trappars_filename = self.get_reference_file(input_model, "trappars")
            self.persat_filename = self.get_reference_file(input_model, "persat")

            # Is any reference file missing?
            missing = False
            missing_reftypes = []
            if self.persat_filename == "N/A":
                missing = True
                missing_reftypes.append("PERSAT")
            if self.trap_density_filename == "N/A":
                missing = True
                missing_reftypes.append("TRAPDENSITY")
            if self.trappars_filename == "N/A":
                missing = True
                missing_reftypes.append("TRAPPARS")
            if missing:
                if len(missing_reftypes) == 1:
                    msg = "Missing reference file type:  " + missing_reftypes[0]
                else:
                    msg = "Missing reference file types: "
                    for name in missing_reftypes:
                        msg += " " + name
                self.log.warning("%s", msg)
                input_model.meta.cal_step.persistence = "SKIPPED"
                return input_model

            # Work on a copy
            result = input_model.copy()

            if self.input_trapsfilled is None:
                traps_filled_model = None
            else:
                traps_filled_model = datamodels.TrapsFilledModel(self.input_trapsfilled)
            trap_density_model = datamodels.TrapDensityModel(self.trap_density_filename)
            trappars_model = datamodels.TrapParsModel(self.trappars_filename)
            persat_model = datamodels.PersistenceSatModel(self.persat_filename)

            pers_a = persistence.DataSet(
                result,
                traps_filled_model,
                self.flag_pers_cutoff,
                self.save_persistence,
                trap_density_model,
                trappars_model,
                persat_model,
            )
            (result, traps_filled, output_pers, skipped) = pers_a.do_all()
            if skipped:
                result.meta.cal_step.persistence = "SKIPPED"
            else:
                result.meta.cal_step.persistence = "COMPLETE"

            if traps_filled_model is not None:  # input traps_filled
                del traps_filled_model
            if traps_filled is not None:  # output traps_filled
                # Save the traps_filled image with suffix 'trapsfilled'.
                self.save_model(traps_filled, "trapsfilled", force=self.save_trapsfilled)
                del traps_filled

            if output_pers is not None:  # output file of persistence
                self.save_model(output_pers, suffix="output_pers")
                output_pers.close()

            # Cleanup
            del trap_density_model
            del trappars_model
            del persat_model

        return result
