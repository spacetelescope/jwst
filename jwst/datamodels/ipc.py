from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base

__all__ = ['IPCModel']

class IPCModel(model_base.DataModel):
    """
    A data model for IPC kernel checking information.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst_lib.models.DataModel`.

    data : numpy array
        The deconvolution kernel (a very small image).
    """
    schema_url = "ipc.schema.yaml"

    def __init__(self, init=None, data=None, **kwargs):
        super(IPCModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data
