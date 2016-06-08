#! /usr/bin/env python

from jwst.stpipe import Step, cmdline
from . import persistence

class PersistenceStep( Step ):
    """
    PersistenceStep: Step for removing persistence from exposures. This
    is currently a no-op step.
    """
    def process( self, input ):

        pers_a = persistence.DataSet( input )
        output_obj = pers_a.do_all()

        if output_obj is not None:
            output_obj.meta.cal_step.persistence = 'SKIPPED' # no-op

        return output_obj


if __name__ == '__main__':
    cmdline.step_script( persistence_step )
