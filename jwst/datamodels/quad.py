from .model_base import JwstDataModel


__all__ = ['QuadModel']


class QuadModel(JwstDataModel):
    """
    A data model for 4D image arrays.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/quad.schema"

    def __init__(self, init=None, **kwargs):
        super(QuadModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
