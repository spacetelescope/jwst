#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from . import persistence

class PersistenceStep(Step):
    """
    PersistenceStep: Correct a science image for persistence.
    This is currently a no-op step.
    """

    spec = """
        # This is the name of the most recent trapsfilled file for the
        # current detector.
        input_traps = string(default="")
    """

    # xxx find a better name than perssat
    # This is currently commented out to prevent CRDS from trying to
    # find these files (because they haven't been delivered yet).
    # xxx xxx reference_file_types = ["trapdensity", "traps", "perssat"]

    def process(self, input):

        # Skip all processing for now ...
        output_obj = datamodels.open(input).copy()
        output_obj.meta.cal_step.persistence = 'SKIPPED'
        self.log.warning('Persistence step is currently a no-op: SKIPPING')

        return output_obj
        # ... end skip all processing for now

        if self.input_traps is not None:
            if self.input_traps == "None" or len(self.input_traps) == 0:
                self.input_traps = None

        output_obj = datamodels.open(input).copy()

        self.trap_density_filename = self.get_reference_file(output_obj,
                                                             'trapdensity')
        self.traps_filename = self.get_reference_file(output_obj, 'traps')
        self.pers_sat_filename = self.get_reference_file(output_obj, 'perssat')

        if self.input_traps is None:
            traps_filled_model = None
        else:
            traps_filled_model = datamodels.CubeModel(self.input_traps)
        trap_density_model = datamodels.ImageModel(self.trap_density_filename)
        traps_model = datamodels.TrapsModel(self.traps_filename)
        pers_sat_model = datamodels.SaturationModel(self.pers_sat_filename)

        pers_a = persistence.DataSet(output_obj, traps_filled_model,
                                     trap_density_model, traps_model,
                                     pers_sat_model)
        (output_obj, traps_filled, skipped) = pers_a.do_all()
        if skipped:
            output_obj.meta.cal_step.persistence = 'SKIPPED'
        else:
            output_obj.meta.cal_step.persistence = 'COMPLETE'

        # Save the traps_filled image, using the input file name but with
        # suffix 'trapsfilled'.
        self.save_model(traps_filled, 'trapsfilled')

        traps_filled_model.close()
        trap_density_model.close()
        traps_model.close()
        pers_sat_model.close()

        return output_obj


if __name__ == '__main__':
    cmdline.step_script(persistence_step)
