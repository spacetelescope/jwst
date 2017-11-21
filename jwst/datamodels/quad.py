from . import model_base


__all__ = ['QuadModel']


class QuadModel(model_base.DataModel):
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

    def __init__(self, init=None, data=None, dq=None, err=None, **kwargs):
        super(QuadModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
