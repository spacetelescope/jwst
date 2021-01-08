"""
JWST-specific Step and Pipeline base classes.
"""
from .step import Step
from .pipeline import Pipeline
from .. import datamodels


class JwstStep(Step):
    @classmethod
    def datamodels_open(cls, init, **kwargs):
        return datamodels.open(init, **kwargs)


class JwstPipeline(Pipeline):
    @classmethod
    def datamodels_open(cls, init, **kwargs):
        return datamodels.open(init, **kwargs)
