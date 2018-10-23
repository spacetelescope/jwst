from .model_base import DataModel

__all__ = ['QuadModel']


class QuadModel(DataModel):
    """
    A data model for 4D image arrays.

    Attributes
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array
    """
    schema_url = "quad.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(QuadModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
