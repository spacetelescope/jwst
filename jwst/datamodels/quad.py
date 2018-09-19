from .model_base import DataModel

__all__ = ['QuadModel']


class QuadModel(DataModel):
    """
    A data model for 4D image arrays.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.  4-D.

    dq : numpy array
        The data quality array.  4-D.

    err : numpy array
        The error array.  4-D
    """
    schema_url = "quad.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(QuadModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
