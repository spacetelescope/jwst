from .reference import ReferenceFileModel


__all__ = ['GainModel']


class GainModel(ReferenceFileModel):
    """
    A data model for 2D gain.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The 2-D gain array
    """
    schema_url = "gain.schema.yaml"
