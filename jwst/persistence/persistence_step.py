#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import persistence

__all__ = ["PersistenceStep"]


class PersistenceStep(Step):
    """
    PersistenceStep: Correct a science image for persistence.
    """

    class_alias = "persistence"

    spec = """
        input_trapsfilled = string(default="") # Name of the most recent trapsfilled file for the current detector
        flag_pers_cutoff = float(default=40.) # Pixels with persistence correction >= this value in DN will be flagged in the DQ
        save_persistence = boolean(default=False) # Save subtracted persistence to an output file with suffix '_output_pers'
        save_trapsfilled = boolean(default=True) # Save updated trapsfilled file with suffix '_trapsfilled'
    """

    reference_file_types = ["trapdensity", "trappars", "persat"]

    def process(self, input):

        if self.input_trapsfilled is not None:
            if (self.input_trapsfilled == "None" or
                    len(self.input_trapsfilled) == 0):
                self.input_trapsfilled = None

        output_obj = datamodels.RampModel(input).copy()

        self.trap_density_filename = self.get_reference_file(output_obj,
                                                             "trapdensity")
        self.trappars_filename = self.get_reference_file(output_obj,
                                                         "trappars")
        self.persat_filename = self.get_reference_file(output_obj, "persat")

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
                    msg += (" " + name)
            self.log.warning("%s", msg)
            output_obj.meta.cal_step.persistence = "SKIPPED"
            return output_obj

        if self.input_trapsfilled is None:
            traps_filled_model = None
        else:
            traps_filled_model = datamodels.TrapsFilledModel(
                self.input_trapsfilled)
        trap_density_model = datamodels.TrapDensityModel(
            self.trap_density_filename)
        trappars_model = datamodels.TrapParsModel(self.trappars_filename)
        persat_model = datamodels.PersistenceSatModel(self.persat_filename)

        pers_a = persistence.DataSet(output_obj, traps_filled_model,
                                     self.flag_pers_cutoff,
                                     self.save_persistence,
                                     trap_density_model, trappars_model,
                                     persat_model)
        (output_obj, traps_filled, output_pers, skipped) = pers_a.do_all()
        if skipped:
            output_obj.meta.cal_step.persistence = 'SKIPPED'
        else:
            output_obj.meta.cal_step.persistence = 'COMPLETE'

        if traps_filled_model is not None:      # input traps_filled
            traps_filled_model.close()
        if traps_filled is not None:            # output traps_filled
            # Save the traps_filled image with suffix 'trapsfilled'.
            self.save_model(
                traps_filled, suffix='trapsfilled', force=self.save_trapsfilled
            )
            traps_filled.close()

        if output_pers is not None:             # output file of persistence
            self.save_model(output_pers, suffix='output_pers')
            output_pers.close()

        # Close reference files.
        trap_density_model.close()
        trappars_model.close()
        persat_model.close()

        return output_obj
