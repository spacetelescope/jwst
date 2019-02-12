"""
This module contains lists of modes grouped in different ways
"""
from ..associations.lib.dms_base import (ACQ_EXP_TYPES, IMAGE2_SCIENCE_EXP_TYPES,
                                         IMAGE2_NONSCIENCE_EXP_TYPES,
                                         SPEC2_SCIENCE_EXP_TYPES)

IMAGING_TYPES = set(tuple(ACQ_EXP_TYPES) + tuple(IMAGE2_SCIENCE_EXP_TYPES)
                    + tuple(IMAGE2_NONSCIENCE_EXP_TYPES) +
                    ('fgs_image', 'fgs_focus'))

SPEC_TYPES = SPEC2_SCIENCE_EXP_TYPES

# FGS guide star exposures
FGS_GUIDE_EXP_TYPES = [
    'fgs_acq1',
    'fgs_acq2',
    'fgs_fineguide',
    'fgs_id-image',
    'fgs_id-stack',
    'fgs_track',
]
