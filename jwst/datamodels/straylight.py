from .reference import ReferenceFileModel


__all__ = ['StrayLightModel']


class StrayLightModel(ReferenceFileModel):
    """
    A data model for 2D straylight mask.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        2-D straylight mask array.
    """
    schema_url = "straylight.schema.yaml"
