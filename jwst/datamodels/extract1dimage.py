from .model_base import JwstDataModel


__all__ = ['Extract1dImageModel']


class Extract1dImageModel(JwstDataModel):
    """
    A data model for the extract_1d reference image array.

    Parameters
    __________
    data : numpy float32 array
         1-D extraction regions array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/extract1dimage.schema"
