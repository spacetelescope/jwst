from . import model_base

__all__ = ['Extract1dImageModel']

class Extract1dImageModel(model_base.DataModel):
    """
    A data model for the extract_1d reference image array.

    Parameters
    ----------
    init: any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data: numpy array
        An array of values that define the extraction regions.
    """
    schema_url = "extract1dimage.schema.yaml"

    def __init__(self, init=None, data=None, **kwargs):
        super(Extract1dImageModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data
