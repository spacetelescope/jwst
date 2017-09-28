#
#  Top level module for 2d extraction.
#
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .nirspec import nrs_extract2d
from .grisms import extract_grism_objects


def extract2d(input_model, which_subarray=None, apply_wavecorr=False, reffile=""):
    nrs_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ', 'NRS_LAMP']
    grism_modes = ['NIS_WFSS', 'NRC_GRISM']
    exp_type = input_model.meta.exposure.type.upper()

    if exp_type in nrs_modes:
        return nrs_extract2d(input_model, which_subarray=None,
                              apply_wavecorr=False, reffile="")

    if exp_type in grism_modes:
        return extract_grism_objects(input_model, grism_objects=[], reffile="")