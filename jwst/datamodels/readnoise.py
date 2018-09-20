from .reference import ReferenceFileModel


__all__ = ['ReadnoiseModel']


class ReadnoiseModel(ReferenceFileModel):
    """
    A data model for 2D readnoise.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        Read noise for all pixels.  2-D.
    """
    schema_url = "readnoise.schema.yaml"
