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

# NIRSPEC lamp mode spec types
NRS_LAMP_MODE_SPEC_TYPES = [
    'brightobj',
    'fixedslit',
    'ifu',
    'msaspec',
]

def is_nrs_lamp(datamodel):
    exp_type = datamodel.meta.exposure.type.lower()
    return exp_type in ['nrs_lamp', 'nrs_autowave']

def is_nrs_linelamp(datamodel):
    lamp_state = datamodel.meta.instrument.lamp_state.lower()
    return lamp_state[0:3] in ['lin', 'ref']

def is_nrs_flatlamp(datamodel):
    lamp_state = datamodel.meta.instrment.lamp_state.lower()
    return lamp_state[0:4] == 'flat'

def is_nrs_slit_linelamp(datamodel):
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    exp_type = datamodel.meta.exposure.type.lower()
    return lamp_mode in ['msaspec', 'fixedslit', 'brightobj'] and is_nrs_linelamp(datamodel) and \
           exp_type in ['nrs_autowave', 'nrs_lamp']

def is_nrs_ifu_linelamp(datamodel):
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == 'ifu' and is_nrs_linelamp(datamodel)

def is_nrs_ifu_flatlamp(datamodel):
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == 'ifu' and is_nrs_flatlamp(datamodel)

def is_nrs_ifu_lamp(datamodel):
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == 'ifu' and is_nrs_lamp(datamodel)

def is_nrs_msaspec_lamp(datamodel):
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == 'msaspec' and is_nrs_lamp(datamodel)

def is_nrs_msaspec_linelamp(datamodel):
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == 'msaspec' and is_nrs_linelamp(datamodel)

def is_nrs_msaspec_flatlamp(datamodel):
    lamp_mode = datamodel.meta.instrument.lamp_mode.lower()
    return lamp_mode == 'msaspec' and is_nrs_flatlamp()

def is_nrs_autoflat(datamodel):
    exp_type = datamodel.meta.exposure.type
    return exp_type.lower() == 'nrs_autoflat'

def is_moving_target(input_models):
    """ Determine if a moving target exposure."""
    model = input_models[0]
    if hasattr(model.meta.target, 'type') and \
        model.meta.target.type is not None and model.meta.target.type.lower() == 'moving':
        return True
    return False
