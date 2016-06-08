#
#  Simple linear pipeline

from jwst.stpipe import LinearPipeline

from jwst.ipc.ipc_step import IPCStep
from jwst.dq_init.dq_init_step import DQInitStep
from jwst.refpix.refpix_step import RefPixStep
from jwst.saturation.saturation_step import SaturationStep
from jwst.dark_current.dark_current_step import DarkCurrentStep
from jwst.linearity.linearity_step import LinearityStep
from jwst.jump.jump_step import JumpStep
from jwst.ramp_fitting.ramp_fit_step import RampFitStep
from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.extract_2d.extract_2d_step import Extract2dStep
from jwst.flatfield.flat_field_step import FlatFieldStep
from jwst.persistence.persistence_step import PersistenceStep
from jwst.straylight.straylight_step import StraylightStep
from jwst.fringe.fringe_step import FringeStep
from jwst.photom.photom_step import PhotomStep

class TestLinearPipeline(LinearPipeline):

    pipeline_steps = [
        ('ipc', IPCStep),
        ('dq_init', DQInitStep),
        ('refpix', RefPixStep),
        ('saturation', SaturationStep),
        ('dark_current', DarkCurrentStep),
        ('linearity', LinearityStep),
        ('jump', JumpStep),
        ('ramp_fit', RampFitStep),
#        ('assign_wcs', AssignWcsStep),
        ('extract_2d', Extract2dStep),
        ('flat_field', FlatFieldStep),
        ('persistence', PersistenceStep),
        ('straylight', StraylightStep),
        ('fringe', FringeStep),
        ('photom', PhotomStep)
        ]


