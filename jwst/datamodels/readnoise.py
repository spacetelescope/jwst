from .reference import ReferenceFileModel


__all__ = ['ReadnoiseModel']


class ReadnoiseModel(ReferenceFileModel):
    """
    A data model for 2D readnoise.

    Parameters
    __________
    data : numpy float32 array
         Read noise
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/readnoise.schema"
