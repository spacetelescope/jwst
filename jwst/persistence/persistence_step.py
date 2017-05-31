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
        # `input_trapsfilled` is the name of the most recent trapsfilled
        # file for the current detector.
        # Pixels that have received a persistence correction greater than
        # or equal to `flag_pers_cutoff` DN will be flagged in the pixeldq
        # extension of the output (rootname_persistence.fits) file.
        # if `save_persistence` is True, the persistence that was
        # subtracted (group by group, integration by integration) will be
        # written to an output file with suffix "_output_pers".
        input_trapsfilled = string(default="")
        flag_pers_cutoff = float(default=40.)
        save_persistence = boolean(default=False)
    """

    # This is currently commented out to prevent CRDS from trying to
    # find these files (because they haven't been delivered yet).
    # xxx xxx reference_file_types = ["trapdensity", "traps", "fullwell"]

    def process(self, input):

        # Skip all processing for now ...
        output_obj = datamodels.open(input).copy()
        output_obj.meta.cal_step.persistence = 'SKIPPED'
        self.log.warning('Persistence step is currently a no-op: SKIPPING')

        return output_obj
        # ... end skip all processing for now

        if self.input_trapsfilled is not None:
            if (self.input_trapsfilled == "None" or
                len(self.input_trapsfilled) == 0):
                self.input_trapsfilled = None

        output_obj = datamodels.open(input).copy()

        self.trap_density_filename = self.get_reference_file(output_obj,
                                                             'trapdensity')
        self.traps_filename = self.get_reference_file(output_obj, 'traps')
        self.full_well_filename = self.get_reference_file(output_obj,
                                                          'fullwell')

        if self.input_trapsfilled is None:
            traps_filled_model = None
        else:
            traps_filled_model = datamodels.CubeModel(self.input_trapsfilled)
        trap_density_model = datamodels.ImageModel(self.trap_density_filename)
        traps_model = datamodels.TrapsModel(self.traps_filename)
        full_well_model = datamodels.SaturationModel(self.full_well_filename)

        pers_a = persistence.DataSet(output_obj, traps_filled_model,
                                     self.flag_pers_cutoff,
                                     self.save_persistence,
                                     trap_density_model, traps_model,
                                     full_well_model)
        (output_obj, traps_filled, output_pers, skipped) = pers_a.do_all()
        if skipped:
            output_obj.meta.cal_step.persistence = 'SKIPPED'
        else:
            output_obj.meta.cal_step.persistence = 'COMPLETE'

        if traps_filled_model is not None:      # input traps_filled
            traps_filled_model.close()
        if traps_filled is not None:            # output traps_filled
            # Save the traps_filled image, using the input file name but
            # with suffix 'trapsfilled'.
            self.save_model(traps_filled, 'trapsfilled')
            traps_filled.close()

        if output_pers is not None:             # output file of persistence
            self.save_model(output_pers, 'output_pers')
            output_pers.close()

        # Close reference files.
        trap_density_model.close()
        traps_model.close()
        full_well_model.close()

        return output_obj


if __name__ == '__main__':
    cmdline.step_script(persistence_step)
