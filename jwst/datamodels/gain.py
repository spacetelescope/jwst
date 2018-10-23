from .reference import ReferenceFileModel


__all__ = ['GainModel']


class GainModel(ReferenceFileModel):
    """
    A data model for 2D gain.

    Attributes
    __________
    data : numpy float32 array
         The gain
    """
    schema_url = "gain.schema.yaml"
