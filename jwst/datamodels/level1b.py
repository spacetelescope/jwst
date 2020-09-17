from .model_base import DataModel

__all__ = ['Level1bModel']


class Level1bModel(DataModel):
    """
    A data model for raw 4D ramps level-1b products.

    Parameters
    __________
    data : numpy uint16 array
         The science data

    zeroframe : numpy uint16 array
         Zeroframe array

    refout : numpy uint16 array
         Reference Output

    group : numpy table
         group parameters table

    int_times : numpy table
         table of times for each integration

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/level1b.schema"

    def __init__(self, init=None, **kwargs):
        super().__init__(init=init, **kwargs)
