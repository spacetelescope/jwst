from .model_base import DataModel

__all__ = ['Extract1dImageModel']

class Extract1dImageModel(DataModel):
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
