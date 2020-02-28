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
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/straylight.schema"
