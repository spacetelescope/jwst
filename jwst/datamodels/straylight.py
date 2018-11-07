from .reference import ReferenceFileModel


__all__ = ['StrayLightModel']


class StrayLightModel(ReferenceFileModel):
    """
    A data model for 2D straylight mask.

    Parameters
    __________
    data : numpy uint8 array
         Straylight mask
    """
    schema_url = "straylight.schema.yaml"
