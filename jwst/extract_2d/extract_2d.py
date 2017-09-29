#
#  Top level module for 2d extraction.
#
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import logging

from .nirspec import nrs_extract2d
from .grisms import extract_grism_objects

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def extract2d(input_model, which_subarray=None, apply_wavecorr=False, reffile=""):
    nrs_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ', 'NRS_LAMP']
    grism_modes = ['NIS_WFSS', 'NRC_GRISM']
    exp_type = input_model.meta.exposure.type.upper()
    log.info('EXP_TYPE is {0}'.format(exp_type))
    
    if exp_type in nrs_modes:
        output_model = nrs_extract2d(input_model, which_subarray=None,
                              apply_wavecorr=False, reffile="")

    elif exp_type in grism_modes:
        output_model = extract_grism_objects(input_model, grism_objects=[], reffile="")

    elif exp_type not in supported_modes:
        log.info("'EXP_TYPE {} not supported for extract 2D".format(exp_type))
        input_model.meta.cal_step.extract_2d = 'SKIPPED'
        return input_model

    else:
        # Set the step status to COMPLETE
        output_model.meta.cal_step.extract_2d = 'COMPLETE'
        del input_model
        return output_model