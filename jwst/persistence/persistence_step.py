#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from . import persistence

class PersistenceStep(Step):
    """
    PersistenceStep: Step for removing persistence from exposures. This
    is currently a no-op step.
    """
    def process(self, input):

        # Skip all processing for now
        output_obj = datamodels.open(input).copy()
        output_obj.meta.cal_step.persistence = 'SKIPPED'
        self.log.warning('Persistence step is currently a no-op: SKIPPING')

        #pers_a = persistence.DataSet(input)
        #output_obj = pers_a.do_all()

        return output_obj


if __name__ == '__main__':
    cmdline.step_script(persistence_step)
