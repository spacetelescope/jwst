#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import guider_cds

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["GuiderCdsStep"]


class GuiderCdsStep (Step):

    """
    This step calculates the countrate for each pixel for FGS modes.
    """

    def process(self, input):
        with datamodels.GuiderRawModel(input) as input_model:
            out_model = guider_cds.guider_cds(input_model)

        out_model.meta.cal_step.guider_cds = 'COMPLETE'

        return out_model
