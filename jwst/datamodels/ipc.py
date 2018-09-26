from .reference import ReferenceFileModel

__all__ = ['IPCModel']

class IPCModel(ReferenceFileModel):
    """
    A data model for IPC kernel checking information.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The deconvolution kernel (a very small image).
    """
    schema_url = "ipc.schema.yaml"
