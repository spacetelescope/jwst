from .reference import ReferenceFileModel

__all__ = ['IPCModel']

class IPCModel(ReferenceFileModel):
    """
    A data model for IPC kernel checking information.

    Parameters
    __________
    data : numpy float32 array
         IPC deconvolution kernel
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/ipc.schema"
