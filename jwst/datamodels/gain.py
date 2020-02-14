from .reference import ReferenceFileModel


__all__ = ['GainModel']


class GainModel(ReferenceFileModel):
    """
    A data model for 2D gain.

    Parameters
    __________
    data : numpy float32 array
         The gain
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/gain.schema"
