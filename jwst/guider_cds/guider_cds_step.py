#! /usr/bin/env python

from crds.core.exceptions import CrdsLookupError

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

    class_alias = "guider_cds"

    def process(self, input):
        with datamodels.GuiderRawModel(input) as input_model:

            # Get the gain reference file
            gain_model = None   # will overwrite if file can be retrieved
            try:
                gain_filename = self.get_reference_file(input_model, 'gain')
                gain_model = datamodels.GainModel(gain_filename)
                self.log.info('Using GAIN reference file: %s', gain_filename)
            except CrdsLookupError:
                self.log.warning('Unable to retrieve GAIN ref file.')

            # Get the readnoise reference file
            readnoise_model = None   # will overwrite if file can be retrieved
            try:
                readnoise_filename = self.get_reference_file(input_model, 'readnoise')
                readnoise_model = datamodels.ReadnoiseModel(readnoise_filename)
                self.log.info('Using READNOISE reference file: %s',
                              readnoise_filename)
            except CrdsLookupError:
                self.log.warning('Unable to retrieve READNOISE ref file.')

            out_model = guider_cds.guider_cds(input_model, gain_model, readnoise_model)

        out_model.meta.cal_step.guider_cds = 'COMPLETE'

        return out_model
