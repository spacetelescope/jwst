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


# JwstPipeline needs to inherit from Pipeline, but also
# be a subclass of JwstStep so that it will pass checks
# when constructing a pipeline using JwstStep class methods.
# It's important that Pipeline occur first so that it
# takes precedence in the resolution order.
class JwstPipeline(Pipeline, JwstStep):
    pass
