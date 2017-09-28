#
#  Top level module for 2d extraction.
#
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .nirspec import nrs_extract2d


def extract2d(input_model, which_subarray=None, apply_wavecorr=False, reffile=""):
    nrs_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ', 'NRS_LAMP']
    exp_type = input_model.meta.exposure.type.upper()

    if exp_type in nrs_modes:
        return nrs_extract2d(input_model, which_subarray=None,
                              apply_wavecorr=False, reffile="")
