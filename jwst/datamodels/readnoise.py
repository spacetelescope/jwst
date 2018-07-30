from . import model_base


__all__ = ['ReadnoiseModel']


class ReadnoiseModel(model_base.DataModel):
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

    def __init__(self, init=None, data=None, **kwargs):
        super(ReadnoiseModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data
