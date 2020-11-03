"""Step parameters model"""
from copy import copy

from .model_base import JwstDataModel


__all__ = ['StepParsModel']


DEFAULT_META = {
    'date': 'SPECIFY DATE',
    'origin': 'STScI',
    'telescope': 'JWST',
    'reftype': 'pars-step',
    'pedigree': 'SPECIFY PEDIGREE',
    'description': 'Parameters for calibration step SPECIFY',
    'author': 'SPECIFY AUTHOR',
    'useafter': 'SPECIFY',
    'instrument': {
        'name': 'SPECIFY',
    },
}


class StepParsModel(JwstDataModel):
    """
    A data model for `Step` parameters.
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/steppars.schema"

    def __init__(self, init=None, **kwargs):
        super().__init__(init=init, **kwargs)
        meta = copy(DEFAULT_META)
        meta.update(self.meta.instance)
        self.meta.instance.update(meta)
