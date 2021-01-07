"""Dummy up the resample/cube_build steps"""

import numpy as np

from ... import datamodels
from ...stpipe import Step


class TempSkipStep(Step):
    """Temporarily disable the parent step

    Notes
    -----
    As a pre-hook, this will set the parent step's `skip` to True so
    the parent step is not executed. The intent is that there is a
    post-hook that will do something instead. The post hook has the
    responsibility of restoring the original skip status which is
    stored in `self.parent._orig_skip`.
    """

    spec = """
    """

    def process(self, *args):
        self.log.warning(
            'Temporarily forcing skipping parent step'
        )
        self.parent._orig_skip = self.parent.skip
        self.parent.skip = True


class MockReflectionStep(Step):
    """Dummy step that simply returns the passe input

    Notes
    -----
    Run as a post hook on the results.
    """

    spec = """
    """

    def process(self, input):
        self.log.warning(
            'Reflecting back the input'
        )

        # Reset the parent's `skip`.
        # If an explicit skip was requested, then skip here also.
        self.parent.skip = getattr(self.parent, '_orig_skip', self.parent.skip)
        if self.parent.skip:
            self.log.warning('Skipping because original parent step is skipped')
            return

        return input


class MockResampleStep(Step):
    """Dummy up results for resample

    Notes
    -----
    Run as a post hook on the results.
    """

    spec = """
    """

    def process(self, input):
        """Create mock resample data
        """
        self.log.warning(
            'Creating mock resampled product until step is available'
        )

        # Reset the parent's `skip`.
        # If an explicit skip was requested, then skip here also.
        self.parent.skip = getattr(self.parent, '_orig_skip', self.parent.skip)
        if self.parent.skip:
            self.log.warning('Skipping because original parent step is skipped')
            return

        # Pick the first model to create the mock data on.
        model = input[0]

        # Mock it
        resampled = datamodels.DrizProductModel()
        resampled.data = model.data
        resampled.wht = model.err
        resampled.con = model.dq
        resampled.meta = model.meta

        # That's all folks
        resampled.meta.cal_step.resample = 'COMPLETE'
        self.log.info('Step completed')
        return resampled


class MockCubeBuildStep(Step):
    """Dummy up results for the cube build step.

    Notes
    -----
    Run as a post hook on the results.
    """

    spec = """
    """

    def process(self, inputs):
        self.log.warning(
            'Creating fake cube build product until step is available'
        )

        # Reset the parent's `skip`.
        # If an explicit skip was requested, then skip here also.
        self.parent.skip = getattr(self.parent, '_orig_skip', self.parent.skip)
        if self.parent.skip:
            self.log.warning('Skipping because original parent step is skipped')
            return

        # Pick the first model to create the mock data on.
        model = input[0]

        # Mock it
        resampled = datamodels.IFUCubeModel()
        naxis = model.data.shape[0]
        array = np.ones((naxis, naxis, naxis))
        resampled.data = array
        resampled.dq = array
        resampled.err = array
        resampled.weightmap = array
        resampled.meta = model.meta

        # That's all folks
        resampled.meta.cal_step.resample = 'COMPLETE'
        self.log.info('Step completed')
        return resampled
