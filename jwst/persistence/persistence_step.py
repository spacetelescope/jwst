#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from . import persistence

class PersistenceStep(Step):
    """
    PersistenceStep: Correct a science image for persistence.
    """

    spec = """
    """

    reference_file_types = ["trapsfilled", "trapdensity", "traps",
                            "gain", "saturation"]

    def process(self, input):

        output_obj = datamodels.open(input).copy()

        self.traps_filled_filename = self.get_reference_file(output_obj,
                                                             'trapsfilled')
        self.trap_density_filename = self.get_reference_file(output_obj,
                                                             'trapdensity')
        self.traps_filename = self.get_reference_file(output_obj, 'traps')
        self.gain_filename = self.get_reference_file(output_obj, 'gain')
        self.sat_filename = self.get_reference_file(output_obj, 'saturation')

        # Is CubeModel appropriate for trap_map_model?
        traps_filled_model = datamodels.CubeModel(self.traps_filled_filename)
        trap_density_model = datamodels.ImageModel(self.trap_density_filename)
        traps_model = datamodels.TrapsModel(self.traps_filename)
        gain_model = datamodels.GainModel(self.gain_filename)
        sat_model = datamodels.SaturationModel(self.sat_filename)

        pers_a = persistence.DataSet(output_obj,
                                     traps_filled_model, trap_density_model,
                                     traps_model, gain_model, sat_model)
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
        gain_model.close()
        sat_model.close()

        return output_obj


if __name__ == '__main__':
    cmdline.step_script(persistence_step)
