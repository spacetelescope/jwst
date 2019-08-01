"""Step parameters model"""
import os
from .model_base import DataModel

__all__ = ['StepParsModel']


class StepParsModel(DataModel):
    """
    A data model for `Step` parameters.
    """
    schema_url = "steppars.schema"
    supported_formats = ['yaml', 'json', 'asdf']

    def __init__(self, init=None, **kwargs):
        super().__init__(init=init, **kwargs)
