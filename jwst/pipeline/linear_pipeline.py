#
#  Simple linear pipeline

from ..stpipe import LinearPipeline

from ..ipc.ipc_step import IPCStep
from ..dq_init.dq_init_step import DQInitStep
from ..refpix.refpix_step import RefPixStep
from ..saturation.saturation_step import SaturationStep
from ..dark_current.dark_current_step import DarkCurrentStep
from ..linearity.linearity_step import LinearityStep
from ..jump.jump_step import JumpStep
from ..ramp_fitting.ramp_fit_step import RampFitStep
from ..assign_wcs.assign_wcs_step import AssignWcsStep
from ..extract_2d.extract_2d_step import Extract2dStep
from ..flatfield.flat_field_step import FlatFieldStep
from ..persistence.persistence_step import PersistenceStep
from ..straylight.straylight_step import StraylightStep
from ..fringe.fringe_step import FringeStep
from ..photom.photom_step import PhotomStep


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
        ('assign_wcs', AssignWcsStep),
        ('extract_2d', Extract2dStep),
        ('flat_field', FlatFieldStep),
        ('persistence', PersistenceStep),
        ('straylight', StraylightStep),
        ('fringe', FringeStep),
        ('photom', PhotomStep)
        ]
