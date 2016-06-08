from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base
from .dynamicdq import dynamic_mask


__all__ = ['FlatModel']


class FlatModel(model_base.DataModel):
    """
    A data model for 2D flat-field images.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst_lib.models.DataModel`.

    data : numpy array
        The science data.  2-D.

    dq : numpy array
        The data quality array.  2-D.

    err : numpy array
        The error array.  2-D.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "flat.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None,
                 dq_def=None, **kwargs):
        super(FlatModel, self).__init__(init=init, **kwargs)

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
