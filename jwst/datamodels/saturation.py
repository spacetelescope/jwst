from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base
from .dynamicdq import dynamic_mask

__all__ = ['SaturationModel']

class SaturationModel(model_base.DataModel):
    """
    A data model for saturation checking information.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst_lib.models.DataModel`.

    data : numpy array
        The science data.

    dq : numpy array
        The data quality array.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "saturation.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, dq_def=None, **kwargs):
        super(SaturationModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if dq_def is not None:
            self.dq_def = dq_def

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
