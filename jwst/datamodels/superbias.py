from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base
from .dynamicdq import dynamic_mask


__all__ = ['SuperBiasModel']


class SuperBiasModel(model_base.DataModel):
    """
    A data model for 2D super-bias images.
    """
    schema_url = "superbias.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None,
                 dq_def=None, **kwargs):
        super(SuperBiasModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        if dq_def is not None:
            self.dq_def = dq_def

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
