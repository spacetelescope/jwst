"""Dummy up the resample/cube_build steps"""

import numpy as np

from ... import datamodels
from ...stpipe import Step


class DummyResampleStep(Step):
    """Dummy up results for resample
    """

    spec = """
    """

    def process(self, input):
        self.log.warning(
            'Creating fake resampled product until step is available'
        )
        resampled = datamodels.DrizProductModel()
        resampled.data = input[0].data
        resampled.wht = input[0].err
        resampled.con = input[0].dq
        resampled.relsens = input[0].relsens
        resampled.meta = input[0].meta

        resampled.meta.cal_step.resample = 'COMPLETE'
        self.log.info('Step completed')
        return resampled


class DummyCubeBuildStep(Step):
    """Dummy up results for
    """

    spec = """
    """

    def process(self, input):
        self.log.warning(
            'Creating fake resampled product until step is available'
        )

        resampled = datamodels.IFUCubeModel()
        naxis = input[0].data.shape[0]
        array = np.ones((naxis, naxis, naxis))
        resampled.data = array
        resampled.dq = array
        resampled.err = array
        resampled.weightmap = array
        resampled.meta = input[0].meta

        resampled.meta.cal_step.resample = 'COMPLETE'
        self.log.info('Step completed')
        return resampled
