from .reference import ReferenceFileModel

__all__ = ['BarshadowModel']


class BarshadowModel(ReferenceFileModel):
    """
    A data model for Bar Shadow correction information.

    Parameters
    __________
    data1x1 : numpy float32 array
         Bar Shadow 1x1 data array

    var1x1 : numpy float32 array
         Bar Shadow 1x1 correction variance

    data1x3 : numpy float32 array
         Bar Shadow 1x3 data array

    var1x3 : numpy float32 array
         Bar Shadow 1x3 correction variance
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/barshadow.schema"
